module SchwarzschildGeo

using DifferentialEquations, Images, ProgressMeter

const black = RGB{N0f8}(0,0,0)
const bh_image_file = joinpath(dirname(@__DIR__),"media","galaxy.jpg")
const DEFAULT_TSPAN = (0.,10000.)
const DEFAULT_R0 = 20.0
const DEFAULT_MAX_R = 10000.

f(r) = 1-1/r
f′(r) = 1/r^2

function nullgeoeom!(dx,x,p,τ)
	@debug "Parameter" τ
	@debug "Coordinate Location" x
	@debug "Parameters" p

	r0 = p[1]
	θ0 = p[2]
	ṫ0 = p[3]
	ϕ̇0 = p[4]

	# x[1] = r, x[2] = ṙ, x[3] = θ, x[4] = θ̇, x[5] = t, x[6] = ϕ
	r = x[1]
	θ = x[3]

	f0 = f(r0)
	fr = f(r)
	fp = f′(r)
	
	ṫ = dx[5] = f0/fr*ṫ0
	ϕ̇ = dx[6] = (r0*sin(θ0)/r/sin(θ))^2*ϕ̇0
	θ̇ = dx[3] = x[4]
	ṙ = dx[1] = x[2]
	#=r̈=#dx[2] = (1/2)*(-fr*fp*ṫ^2 + fp/fr*ṙ^2 + 2*r*fr*(θ̇^2+sin(θ)^2*ϕ̇^2))
	#=θ̈=#dx[4] = (1/2)*sin(2*θ)*ϕ̇^2 - 2*ṙ/r*θ̇
	
end

# x⃗  = (r, θ, t, ϕ)
# r⃗ᵥ = (ϕv, θv)
function first_order_initial_from_view(x0, ang_view)
	ϕv, θv = ang_view
	r₀, θ₀ = x0
	sqf₀ = sqrt(f(r₀))
	return (-sqf₀*sin(θv)*cos(ϕv), cos(θv)/r₀, 1/sqf₀, sin(θv)*sin(ϕv)/r₀/sin(θ₀))
end


# x⃗  = (r, θ, t, ϕ)
# r⃗ᵥ = (ϕv, ϕv)
function sch_ode_make(x0, ang_view; tspan=DEFAULT_TSPAN)

	# r0 = p[1]
	# θ0 = p[2]
	# ṫ0 = p[3]
	# ϕ̇0 = p[4]
	ẋ₀ = first_order_initial_from_view(x0, ang_view)
	p = [x0[1], x0[2], ẋ₀[3], ẋ₀[4]]

	# x[1] = r, x[2] = ṙ, x[3] = θ, x[4] = θ̇, x[5] = t, x[6] = ϕ
	x0dx0 = [x0[1], ẋ₀[1], x0[2], ẋ₀[2], x0[3], x0[4]]

	return ODEProblem(nullgeoeom!, x0dx0, tspan, p)
end


function sch_ode_solve(ϕv, θv; r0=DEFAULT_R0, θ0 = pi/2, t0=0.0, ϕ0=0.0, 
											 tspan=DEFAULT_TSPAN, outerR = DEFAULT_MAX_R)
	x0 = (r0, θ0, t0, ϕ0)
	
	ode = sch_ode_make(x0, (ϕv, θv); tspan)
	
	condition(x,t,integrator) = x[1] < 1 || x[1] > outerR
	cb = DiscreteCallback(condition, terminate!)

	sol = solve(ode, callback=cb)
	return sol
	
end

function schview(ϕv, θv; r0=DEFAULT_R0, θ0 = pi/2, ϕ0=0.0, tspan=DEFAULT_TSPAN)
	sol = sch_ode_solve(ϕv, θv; r0=r0, θ0=θ0, ϕ0=ϕ0, tspan=tspan)
	
	final = sol[end]

	if final[1]<=1
		return :black, :black
	end
	
	ϕout, θout = final[[6,3]]

	return ϕout, θout 
end

function angletopix(ϕout, θout, imgdim)
	y = mod(θout, 2*pi)
	y = y < pi ? y : 2*pi-y
	y = y/pi*imgdim[1]
	yp = y |> ceil |> Int

	xp = mod(ϕout, 2*pi)/(2*pi)*imgdim[2] |> ceil |> Int

	if yp <= 0 
		yp=1
	elseif yp > imgdim[1]
		yp=imgdim[1]
	end

	if xp == 0 
		xp+=1
	end

	yp, xp
end

function schpic(imgfile=bh_image_file; r0 = DEFAULT_R0, tspan=DEFAULT_TSPAN, 
		ϕcam=pi, θcam=0.0, ϕ0=0.0)

	img = load(imgfile)
	imgout = similar(img)

	h,w = size(img)

	@showprogress	for j in 1:h, i in 2:w-1
		ϕin = i/w*2*pi 
		θin = j/h*pi   
		
		if θin ≈ 0.0
			θin += j/h*pi
		elseif θin ≈ pi
			θin -= j/h*pi
		end

		ϕout, θout = schview(ϕin - ϕcam, θin - θcam; r0=r0, tspan=tspan, ϕ0=ϕ0)

		if ϕout == :black
			imgout[j,i] = black
			continue
		end

		viewpixind = angletopix(ϕout, θout, (h,w))
		imgout[j,i] = img[viewpixind...]
	end

	imgout
end

export schview, schpic, sch_ode_make, sch_ode_solve, bh_image_file

end # module
