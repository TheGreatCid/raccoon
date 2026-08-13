#include "InversionSafeDamper.h"

#include "FEProblemBase.h"
#include "DisplacedProblem.h"
#include "MooseMesh.h"
#include "Assembly.h"

#include "libmesh/elem.h"
#include "libmesh/node.h"

registerMooseObject("raccoonApp", InversionSafeDamper);

InputParameters
InversionSafeDamper::validParams()
{
  InputParameters params = GeneralDamper::validParams();
  params.addClassDescription(
      "Damps the Newton step so that no element inverts on the displaced mesh. The step-to-inversion "
      "is found algebraically from det(J_map(lambda)) = det(J0 + lambda*dJ), so it caps an inverting "
      "step instead of failing on it (as ElementJacobianDamper does).");
  params.addRequiredParam<std::vector<VariableName>>("displacements",
                                                     "The displacement variables");
  params.addParam<Real>("safety_factor",
                        0.9,
                        "Allow at most this fraction of the step to the first element inversion.");
  params.addParam<Real>(
      "min_relative_jacobian",
      0.2,
      "Limit the step before any element's FE-map Jacobian falls below this fraction of its "
      "current (lambda=0) value.");
  params.addParam<unsigned int>(
      "n_samples", 32, "Number of samples along the step used to bracket the first inversion.");
  params.addParam<unsigned int>(
      "quadrature_order",
      0,
      "Order of the Gauss rule used to sample the element map (0 = the FE type's default). Set it "
      "to match the [Quadrature] order so the same points the assembly integrates on are guarded.");
  return params;
}

InversionSafeDamper::InversionSafeDamper(const InputParameters & parameters)
  : GeneralDamper(parameters),
    _tid(parameters.get<THREAD_ID>("_tid")),
    _fe_problem(*parameters.getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")),
    _safety_factor(getParam<Real>("safety_factor")),
    _min_relative_jacobian(getParam<Real>("min_relative_jacobian")),
    _n_samples(getParam<unsigned int>("n_samples"))
{
  // Work on the undisplaced (reference) mesh: it carries the displacement dofs and its node
  // coordinates are the reference coordinates x_ref. The deformed configuration is reconstructed
  // from the solution/update vectors, so we never rely on the displaced mesh's sync state.
  _mesh = &_fe_problem.mesh();
  _dim = _mesh->dimension();

  const auto & disp_names = getParam<std::vector<VariableName>>("displacements");
  _ndisp = disp_names.size();
  for (const auto & name : disp_names)
    _disp_var.push_back(&_sys.getFieldVariable<Real>(_tid, name));

  // A standalone FE + quadrature to obtain the (geometry-independent) master-element shape
  // gradients. Use the requested order, or the FE type's default. (The assembly's own quadrature
  // is not yet available at construction, so we do not query it here.)
  libMesh::FEType fe_type = _disp_var[0]->feType();
  const unsigned int req_order = getParam<unsigned int>("quadrature_order");
  const libMesh::Order qorder =
      req_order > 0 ? static_cast<libMesh::Order>(req_order) : fe_type.default_quadrature_order();
  _fe = libMesh::FEBase::build(_dim, fe_type);
  _qrule = std::make_unique<libMesh::QGauss>(_dim, qorder);
  _fe->attach_quadrature_rule(_qrule.get());
  _dphidxi = &_fe->get_dphidxi();
  _dphideta = (_dim >= 2) ? &_fe->get_dphideta() : nullptr;
  _dphidzeta = (_dim >= 3) ? &_fe->get_dphidzeta() : nullptr;
}

Real
InversionSafeDamper::computeDamping(const NumericVector<Number> & solution,
                                    const NumericVector<Number> & update)
{
  Real min_lambda = 1.0;

  PARALLEL_TRY
  {
    auto range = _mesh->getMesh().active_local_element_ptr_range();
    if (range.begin() != range.end())
    {
      // Reinit once on a (valid) current element to populate the master-element shape gradients,
      // which are the same for all elements of this type. We never reinit a folded configuration.
      _fe->reinit(*range.begin());
      const unsigned int nqp = _qrule->n_points();
      const unsigned int nsh = _dphidxi->size();

      std::vector<std::vector<RealVectorValue>> G(nsh, std::vector<RealVectorValue>(nqp));
      for (unsigned int i = 0; i < nsh; ++i)
        for (unsigned int q = 0; q < nqp; ++q)
        {
          G[i][q](0) = (*_dphidxi)[i][q];
          if (_dphideta)
            G[i][q](1) = (*_dphideta)[i][q];
          if (_dphidzeta)
            G[i][q](2) = (*_dphidzeta)[i][q];
        }

      for (auto & elem : range)
      {
        const unsigned int nn = elem->n_nodes();

        // Reconstruct the OLD (current) deformed node coordinates and the Newton search direction
        // purely from the vectors. MOOSE calls this with soln = the full undamped new solution
        // (= old - update), so the old displacement is (solution + update); the deformed step at
        // scale lambda is x_ref + (solution + update) - lambda*update, i.e. x_old - lambda*du.
        std::vector<Point> x(nn), du(nn);
        for (unsigned int i = 0; i < nn; ++i)
        {
          const Node & node = elem->node_ref(i);
          x[i] = node; // reference coordinate x_ref
          for (unsigned int j = 0; j < _ndisp; ++j)
          {
            const dof_id_type dof = node.dof_number(_sys.number(), _disp_var[j]->number(), 0);
            const Real upd = update(dof);
            x[i](j) += solution(dof) + upd; // x_old = x_ref + old displacement
            du[i](j) = upd;                 // search direction (subtracted below)
          }
        }

        for (unsigned int q = 0; q < nqp; ++q)
        {
          // J(lambda) = J0 - lambda*dJ; fill any unused dimension with identity.
          RealTensorValue J0, dJ;
          for (unsigned int d = _dim; d < 3; ++d)
            J0(d, d) = 1.0;
          for (unsigned int i = 0; i < nn && i < nsh; ++i)
            for (unsigned int a = 0; a < _dim; ++a)
              for (unsigned int b = 0; b < _dim; ++b)
              {
                J0(a, b) += x[i](a) * G[i][q](b);
                dJ(a, b) += du[i](a) * G[i][q](b);
              }

          const Real det0 = J0.det();
          if (det0 <= 0.0)
          {
            // Current configuration already (nearly) inverted -- creep forward.
            min_lambda = std::min(min_lambda, Real(1e-3));
            continue;
          }
          const Real floor = _min_relative_jacobian * det0;

          // Find the first lambda in (0,1] where det drops to the floor, then bisect.
          Real prev = 0.0;
          for (unsigned int s = 1; s <= _n_samples; ++s)
          {
            const Real l = Real(s) / Real(_n_samples);
            if ((J0 - l * dJ).det() < floor)
            {
              Real lo = prev, hi = l;
              for (unsigned int it = 0; it < 40; ++it)
              {
                const Real mid = 0.5 * (lo + hi);
                if ((J0 - mid * dJ).det() < floor)
                  hi = mid;
                else
                  lo = mid;
              }
              min_lambda = std::min(min_lambda, lo);
              break;
            }
            prev = l;
          }
        }
      }
    }
  }
  PARALLEL_CATCH;

  _communicator.min(min_lambda);

  // Take the full Newton step unless a fold was actually bracketed within it (min_lambda < 1).
  // Only then back off to `safety_factor` of the distance-to-fold. Applying safety_factor to an
  // un-reduced min_lambda would needlessly shrink every step and stall convergence.
  Real damping = 1.0;
  if (min_lambda < 1.0)
  {
    damping = std::max(std::min(_safety_factor * min_lambda, Real(1.0)), Real(1e-6));
    _console << "InversionSafeDamper '" << name()
             << "': limiting step to prevent inversion, damping = " << damping
             << " (lambda_to_fold = " << min_lambda << ")" << std::endl;
  }
  return damping;
}
