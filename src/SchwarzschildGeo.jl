module SchwarzschildGeo

using DifferentialEquations

function radial_U_diff!(du,u,p,ϕ)
	du[1] = u[2]
	du[2] = 3*(u[1])^2 - u[1]
end

function initial_U_cond(θin, rin=20)
	[
	 1/rin,
	 1/(rin*tan(θin))
	 ]
end

sch_f(r) = 1-2/r

L(θv, r0=20) = r0*sin(θv)
vr0(θv, r0=20, E=1) = -√(E^2 - sch_f(r0)*sin(θv))

function radial_R_diff!(dr,r,p,τ)
	ℓ = L(p[1],p[2])
	# r[2] ≡ ṙ
	# r[1] ≡ r
	dr[1] = r[2]
	dr[2] = ℓ^2*(r[1]^-3. - 3*r[1]^-4.)
	dr[3] = ℓ/(r[1])^2.0
end

function sch_ode_prob(θv; r0=20., tspan=(0.,10.))
	ODEProblem(radial_R_diff!, [r0, vr0(θv, r0), 0.0], tspan, [θv,r0])
end 

export radial_U_diff!, initial_U_cond, sch_ode_prob

end # module
