module SchwarzschildGeo

using DifferentialEquations, Images

const black = RGB{N0f8}(0,0,0)
const bh_image_file = joinpath(dirname(@__DIR__),"media","galaxy.jpg")

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

L(ϕv, r0=20) = sin(ϕv)
vr0(ϕv, r0=20, E=1) = -cos(ϕv)#-√(E^2*r0 - sch_f(r0)*sin(ϕv))

function radial_R_diff!(dr,r,p,τ)
	ℓ = L(p[1],p[2])
	# r[2] ≡ ṙ
	# r[1] ≡ r
	dr[1] = r[2]
	dr[2] = ℓ^2*(r[1]^-3. - 3*r[1]^-4.)
	dr[3] = ℓ/(r[1])^2.0
end

function sch_ode_prob(ϕv; r0=20., tspan=(0.,10.))
	ODEProblem(radial_R_diff!, [r0, vr0(ϕv, r0), ϕv], tspan, [ϕv,r0])
end 

function sch_ode_solve(ϕv; r0=20., tspan=(0.,10.))
	condition(r,t,integrator) = r[1]>500 || r[1]<2
	cb=DiscreteCallback(condition, terminate!)
	solve(sch_ode_prob(ϕv; r0=r0, tspan=tspan), callback=cb)
end 

function tophi(ϕp,θp)
	sqt = sqrt(cos(θp)^2 + sin(ϕp)^4)

	cψ = sin(ϕp)^2/sqt
	sψ = cos(θp)/sqt

	return atan(cψ*sin(θp)*sin(ϕp) + sψ*cos(θp), sin(ϕp)*cos(ϕp)) + π, cψ, sψ
	
	ϕp, 0, 0
end


function fromphi(ϕ, cψ, sψ)
	sψ *= -1

	ϕp = atan(cψ*tan(ϕ),1)
	θp = atan(sqrt(cos(ϕ)^2+cψ^2*sψ^2)/sψ*sin(ϕ)) + pi/2

	return (ϕp, θp)
end

function angletopix(ϕp, θp, imgdim)
	x = mod(ϕp, 2*pi)/(2*pi)*imgdim[2] |> ceil |> Int

	y = (pi - mod(θp, pi))/pi*imgdim[1]  |> ceil |> Int

	y,x
end

function schview(ϕin, θin; r0=20., tspan=(0.,100.))
	# ϕ, cψ, sψ = tophi(ϕin, θin)
	sol = sch_ode_solve(ϕin; r0=r0, tspan=tspan)
	
	final = sol[end]

	if final[1]<=2
		return :black, :black
	end
	
	ϕsol = final[3]
	ϕout, θout = ϕsol, θin
	# ϕout, θout = fromphi(ϕsol, cψ, sψ)

	return ϕout, θout 
end

function schpic(imgfile=bh_image_file; r0 = 20.0, tspan=(0,100))
	img = load(imgfile)
	imgout = similar(img)

	h,w = size(img)
	
	for j in 1:h, i in 1:w
		# img[j,i]
		ϕin = i/w*2*pi
		θin = (j-h)/h*pi
		
		ϕout, θout = schview(ϕin, θin; r0=r0, tspan=tspan)

		if ϕout == :black
			imgout[j,i] = black
			continue
		end

		viewpixind = angletopix(ϕout+pi, θout, (h,w))
		
		imgout[j,i] = img[viewpixind...]
	end

	imgout
end

export radial_U_diff!, initial_U_cond, sch_ode_prob, schview, schpic

end # module
