# Independent formulation review

Reviewed `solve.py` on 2026-10-13. This review concerns mathematical implementation and interpretation; it is not experimental validation.

## Formulation checks

- The calculation is a genuine three-dimensional tetrahedral finite-element calculation. All three velocity components, pressure and temperature are nodal P1 fields. SUPG momentum stabilization and PSPG pressure stabilization make the equal-order choice intentional.
- Momentum pressure coupling has the expected weak sign, `-p div(v)`. The PSPG pressure-gradient term and right-hand-side body-force term have matching positive signs.
- Momentum and PSPG pseudo-time terms both include the old-velocity contribution, so the extra pseudo-time terms vanish at a fixed point. The separately recomputed final steady nonlinear residual is a useful verification beyond the update criterion.
- The grad-div block indices correctly represent `gamma div(u) div(v)`. Gamma = 0.5 m²/s is a numerical penalty, much larger than the assumed effective viscosity. Its effect should be included in a discretization sensitivity study; it does not make the P1 velocity locally divergence-free.
- The open-boundary backflow matrix uses a nonnegative coefficient `0.5 max(-u.n,0)` and is the usual kinetic-energy stabilization for a natural traction boundary.
- The final heat transport uses the conservative weak form `-theta u·grad(v) + boundary (u.n) theta v`, with theta measured above ambient. Its boundary term is essential to physical energy conservation. The heat graph diffusion makes off-diagonal entries nonpositive by adding a symmetric edge matrix with zero row sums. With nonnegative source loads and the ambient-temperature Dirichlet nodes, it is a sound bound-preserving algebraic diffusion construction. It contributes mesh-dependent artificial diffusion, so peak temperature and plume spreading require mesh sensitivity. This final heat formulation is not thermal SUPG.
- The natural pressure boundary supplies the pressure reference; an additional pressure pin is not required. The plotted pressure is the modeled room-side finite-element traction-reference pressure, not the closed-door leakage pressure drop or the fan static-pressure capability.
- Prescribing the hood velocity with a normalized smooth profile gives the specified net hood flow. The area-average speed is Q / 0.54 m², while the nodal profile peak is necessarily higher.
- The source is normalized in the quadrature rule to 1200 W total. A separate review confirmed nonnegative assembled P1 load entries. Negative weights in tetrahedral order-3 quadrature therefore did not produce negative thermal loads in these meshes.

## Important reporting checks

1. Report nonlinear residual, thermal residual, iteration convergence, net boundary flux imbalance and the local divergence norm separately. A machine-small net flux error does not establish local incompressibility or mesh accuracy.
2. The previous nonconservative heat-advection formulation produced significant energy errors when velocity was only weakly divergence-free. It was replaced by conservative thermal transport with the correct boundary flux. The final direct heat residual must compare hood advective output plus wall heat loss plus door reaction directly against 1200 W, without subtracting a divergence correction. The separately exported divergence diagnostic describes what the earlier nonconservative form would have introduced; it is not part of the final heat balance. Archived earlier results are not final validated outputs.
3. Temperature is fixed to ambient on every open-door node. This is physically appropriate for makeup inflow; if any significant doorway outflow occurs, the same condition acts as an ambient-temperature sink on that outflow. Check and disclose gross outward flux, which is already exported. In the older both-open run outward doorway flow was negligible, so this caveat did not materially affect that case.
4. The uniform assumed eddy viscosity is a phenomenological closure, not a solved RANS turbulence model. The steady prediction is a conditional scenario comparison, not a prediction of instantaneous fluctuations, smoke capture efficiency, oil droplets or combustion safety.
5. Finite-element streamline pictures are derived from the computed velocity field; their visual smoothness does not indicate mesh convergence.

## Geometry and operating assumptions

- Room 3.3 × 1.65 × 2.2 m; worktop 0.9 m; hood underside 1.6 m; upper envelope 0.6 m tall. The inferred left cabinet width is 0.6 m. Burner centers follow the drawing offsets at x = 0.85 and 1.25 m, y = 0.30 m.
- The door width and height are user-confirmed; exact locations in the plan, hood effective inlet depth/area, appliance envelope and burner sensible heat remain modeling assumptions. Door leaf panels and outdoor wind are omitted.
- Maximum free delivery 0.5 m³/s and maximum static pressure 1300 Pa are distinct points, not simultaneous operating conditions. The quadratic interpolation and orifice coefficient are assumptions, and exhaust-duct/filter losses are omitted.
- For both closed doors, the 5 mm assumed bottom gaps have total physical area 0.009 m². The network predicts about 0.226 m³/s (13.55 m³/min) under the assumed curve. This result is highly dependent on the unknown gap area, discharge coefficient, actual fan curve and duct resistance.
- Replacing each 5 mm gap by a 0.10 m inlet band preserves net volume flow but expands area by 20. At similar profile, velocity falls by 20 and momentum flux by 20. Consequently the closed-case near-gap speed, floor jet and detailed circulation cannot be interpreted quantitatively as a resolved 5 mm-gap prediction. The network flow estimate and broad temperature trends are the more useful outputs.
- If local temperature rise approaches 45 K, beta × deltaT is about 0.15. The constant-property Boussinesq approximation then has a material local-density approximation near the heat source; bulk room temperature estimates are less affected.
- Do not infer smoke capture from all heat reaching the hood or from a single-exit steady tracer balance. With only inward makeup openings, every steady conserved tracer is forced to leave through the hood. Exposure or finite-time capture needs an explicitly defined tracer/source/transient model.

The necessary conservative heat correction was reviewed and its weak-form signs and flux structure are correct. No remaining blocking sign or algebraic assembly error was found in the current formulation. The remaining concerns are numerical resolution, closure assumptions and the explicitly homogenized closed-door case.
