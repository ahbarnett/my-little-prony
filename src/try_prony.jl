# try Prony's method for real combo of real decaying exponentials.
# Barnett 7/21/25
# Ref: T. Sauer, Oberwolfach Snapshot 2018-004-EN
# Prony's method: an old trick for new problems
# https://publications.mfo.de/handle/mfo/1338

using GLMakie, LinearAlgebra, PolynomialRoots, Printf, Random
verb = 1

# make data = y values on grid tt
h = 0.1    # spacing of data in t
n = 10       # samples
tt = (0:n-1)*h    # reg time sample grid
re = [1.0, 3.0, 5.0]   # exp decay rates ("e" denotes exact)
R = length(re)   # num terms (rates)
rng = Xoshiro(43);   # fix seed
ae = randn(rng, R)  # reproducible real amplitudes
println("Try Prony: true decay rates ", re, "\n  ampls ",ae)
f(t) = sum(a*exp(-r*t) for (a,r) in zip(ae,re))   # function
y = f.(tt)     # meas data
noi = 1e-8
y .+= noi*randn(rng, n)   # add noise

# Prony step 1: extract decay rates
K = 6  # num coefs in poly p, ie deg(p)-1, must be >R
M = n-K+1   # resulting num of annihilation conditions
if M<=R @warn "num conditions does not exceed R!"; end
H = [y[i+j-1] for i=1:M, j=1:K]   # Hankel matrix
U,S,Vt = svd(H; full=true)        # also get full nullspace in Vt
@printf "\nn=%d, K=%d, M=%d. \t Hankel sing vals:\n" n K M
display(S)
trunc = 1e2*noi    # sing val cutoff
ii = findall(S .< trunc)
if isempty(ii) @warn @sprintf("no sing vals below trunc=%.3g",trunc); end
i = ii[1]  #min(M,K)    # col index of V to use to get a R sing vec
# unclear whether to use first sing val below trunc, or dip into true nullspace?
# i = K   # ?
v = Vt[:,i]    # R sing vec
z = roots(v)      # the nonlin step
zj = z[findall(real(z)>0 && real(z)<1 && abs(imag(z)) < 0.1 for z in z)]  # find roots nr [0,1]
zj = real.(zj)       # enforce reality
rj = -log.(zj)/h    # extracted rates
nrj = length(rj)
println("found ",nrj," +ve Re roots with fitted rates:\n  ",rj)
if nrj==R
	@printf "\tmax root error %.3g\n" norm(rj - re, Inf)
end
# Prony step 2: lin solve for coeffs given known rates
V = [exp(-r*t) for t in tt, r in rj]   # Vandermonde matrix
aj = V\y
println("fitted amplitudes:\n  ",aj)
yfit = V*aj
@printf "\tcond(V)=%.3g, \tmax resid %.3g\n" cond(V) norm(y - yfit, Inf)

if verb>0
	fig = Figure(size=(500,1000))
	ax1 = Axis(fig[1,1]; xlabel="t", ylabel="y")
	scatterlines!(tt,y; label="data")
	for k=1:R
		lines!(tt,ae[k]*exp.(-re[k]*tt); label=@sprintf("component %d",k))
	end
	scatterlines!(tt,yfit; label="fit")
	axislegend(position=:rb)
	ax2 = Axis(fig[2,1]; xlabel="Re z", ylabel="Im z")
	scatter!(exp.(-h*re),zeros(R), marker=:xcross, markersize=15, label="true roots (from rates)")
	scatter!(real.(z), imag.(z), label="poly roots")
	axislegend()
	display(fig)
end
