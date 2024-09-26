# This is the DBR system. 03/20/2023
# The structure is like:
# | 10 air | Nxleft DBR_left | Nx1 cav 1 | Nxmid DBR_mid | Nx2 cav 2| Nxright DBR_right | 10 air |

using LinearAlgebra
using StaticArrays
using NPZ
using MAT

const x0 = 0.47 # detector location, normalized by L1
const Dmax = 0.2
const s = 0.0041
#===============simulation setup=================#
const resolution =4
const dx = 1/resolution
const Courant = 1
const dt = dx*Courant
const Simtime = 100000*(820)
const Nt = round(Int,Simtime/dt)

#===============physical parameters==============#
#---------------cavity parameters----------------#
const L1 = 2050                     # length of cavity 1
const L2 = 1340                   # length of cavity 2
const Nx1 = round(Int,L1*resolution)
const Nx2 = round(Int,L2*resolution)

const n1 = 3.4                   # AlGaAs, refractive index of nonlinear cavity
const n2 = 3.67                  # GaAs, refractive index of passive cavity
const sigma1 = 0                 # conductivity of cavity 1, MEEP unit
const sigma2 = s/L1
println(sigma2)
#---------------DBR parameters-------------------#
const N_left = 3                       # DBR periods number, left
const N_mid = 6                        # DBR periods number, mid
const N_right = 3                      # DBR periods number, right

const n_AlAs = 3.0                     # AlAs, refractive index of DBR
const n_glass = 1.5                    # glass,refractive index of DBR substrate

const dm1 = 80                         # mid DBR, width of AlAs
const dm2 = 130                         # mid DBR, width of glass
const dl1 = 70                         # left DBR, width of AlAs
const dl2 = 140                         # left DBR, width of glass
const dr1 = 70                         # right DBR, width of AlAs
const dr2 = 140                         # right DBR, width of glass
const Nxdm1 = round(Int,dm1*resolution)
const Nxdm2 = round(Int,dm2*resolution)
const Nxdl1 = round(Int,dl1*resolution)
const Nxdl2 = round(Int,dl2*resolution)
const Nxdr1 = round(Int,dr1*resolution)
const Nxdr2 = round(Int,dr2*resolution)
const Nxair = 10;

const nleft_period = [n_AlAs*ones(Nxdl1);n_glass*ones(Nxdl2)]
const nmid_period = [n_glass*ones(Nxdm2);n_AlAs*ones(Nxdm1)]
const nright_period = [n_glass*ones(Nxdr2);n_AlAs*ones(Nxdr1)]

const Nxleft = N_left*(Nxdl1+Nxdl2)
const Nxmid = N_mid*(Nxdm1+Nxdm2)+Nxdm2
const Nxright = N_right*(Nxdr1+Nxdr2)
#---------------construct index arrays------------#
const Nx = Nxair + Nxleft + Nx1 + Nxmid + Nx2 + Nxright + Nxair
const start = 1+Nxair+Nxleft
const end_c = Nxair+Nxleft+Nx1

n0 = ones(Nxair)
for i=1:N_left;append!(n0,nleft_period);end
append!(n0,n1*ones(Nx1))
for i=1:N_mid;append!(n0,nmid_period);end
append!(n0,n_glass*ones(Nxdm2))
append!(n0,n2*ones(Nx2))
for i=1:N_right;append!(n0,nright_period);end
append!(n0,ones(Nxair))

matwrite("index.mat", Dict("n" => n0, "dx" => dx)); 

const n = n0# refractive array
const sigma = [zeros(Nxair);zeros(Nxleft);sigma1*ones(Nx1);zeros(Nxmid);sigma2*ones(Nx2);zeros(Nxright);zeros(Nxair)]
const spd = 1 ./(n.^2)
const fac1 = 1 .+ (sigma*dt)/2
const fac2 = 1 .- (sigma*dt)/2
#---------------nonlinear gain setting------------#
const omega_a = 2*pi/820
const gamma_perp = 10/(3e5)
const gamma_parallel = 0.001/(3e5)
const D0_array = [zeros(Nxair);zeros(Nxleft);0.5*Dmax*(1 .-cos.(((1:Nx1).-0.5)*2pi/Nx1));zeros(Nxmid);zeros(Nx2);zeros(Nxright);zeros(Nxair)]* gamma_parallel# pumping strength profile

const n_D = round(Int, x0*Nx1)
global detection = round(Int, Nxair+Nxleft+n_D)
# initializing the field

global E = zeros(Nx)
global B = zeros(Nx)
global D = zeros(Nx)
global P_r = zeros(Nx)
global P_i = zeros(Nx)

######### initialize the simulator with a gaussian pulse inside the gain cavity ########
# E=0.001*exp.(-((1:Nx).-detection).^2/200^2)
# D = D0_array/gamma_parallel
# matwrite("ini_E.mat", Dict("data" => E))
# matwrite("ini_D.mat", Dict("data" => D))
# matwrite("dx.mat", Dict("dx" => dx))

######### initialize the simulator with data from ########
file = matopen("xE.mat") 
E = read(file,"data")
close(file)
file = matopen("xB.mat") 
B = read(file,"data")
close(file)
file = matopen("xD.mat") 
D = read(file,"data")
close(file)
file = matopen("xP_r.mat") 
P_r = read(file,"data")
close(file)
file = matopen("xP_i.mat") 
P_i = read(file,"data")
close(file)


# file = matopen("initials.mat") 
# E[1:Nx] = read(file,"ini_E")
# B[1:Nx] = read(file,"ini_B")                   
# D[(1+10+Nxleft):(10+Nxleft+Nx1)] = read(file,"ini_D")
# P_r[(1+10+Nxleft):(10+Nxleft+Nx1)] = read(file,"ini_Pr")
# P_i[(1+10+Nxleft):(10+Nxleft+Nx1)] = read(file,"ini_Pi")


global changeling = @MMatrix zeros(3,3)
changeling[1,1] = -gamma_parallel
changeling[1,2] = 0
changeling[2,1] = 0
changeling[2,2] = -gamma_perp
changeling[2,3] = omega_a
changeling[3,2] = -omega_a
changeling[3,3] = -gamma_perp

# detection data setting
const Sampling = 20
const Nt_sample = round(Int,Nt/Sampling)
global E_data = zeros(Nt_sample+1)
global D_data = zeros(Nt_sample+1)
E_data[Nt_sample+1] = dt*Sampling
D_data[Nt_sample+1] = dt*Sampling




function update_E(E::Vector{Float64},B::Vector{Float64},P_r::Vector{Float64},P_i::Vector{Float64})
    E_old = E[2];
      @inbounds for xii = (2:Nx)
        E[xii] = (fac2[xii]*E[xii] + spd[xii]*(Courant*(B[xii] - B[xii-1]) - 2*dt*(omega_a*P_i[xii] - gamma_perp*P_r[xii])))/fac1[xii] 
      end
    # absorbing left BC
    E[1] = E_old + (Courant-1)/(Courant+1)*(E[2]-E[1])
end
  
function update_B(E::Vector{Float64},B::Vector{Float64})
    B_old = B[Nx-1];
    @inbounds for xii = (1:Nx-1)
      B[xii] += Courant*(E[xii+1] - E[xii])
    end

    # absorbing right BC
    B[Nx] = B_old+(Courant-1)/(Courant+1)*(B[Nx-1]-B[Nx]);
end

function update_uvec(E::Vector{Float64},D::Vector{Float64},P_r::Vector{Float64},P_i::Vector{Float64},changeling::MMatrix{3, 3, Float64, 9})
    @inbounds for xii = (start:end_c)
        changeling[1,3] = gamma_parallel*real(E[xii])
        changeling[3,1] = - gamma_perp * real(E[xii])
        Left_matrix = I - dt*changeling/2
        Right_matrix = I + dt*changeling/2
        uvec0 = @MArray [D[xii];P_r[xii];P_i[xii]]
        B_vec = @MArray [D0_array[xii];0;0]
        uvec1 = Left_matrix\(Right_matrix*uvec0 + dt*B_vec)
        D[xii]=uvec1[1]
        P_r[xii]=uvec1[2]
        P_i[xii]=uvec1[3]
    end
end

function run(E::Vector{Float64},B::Vector{Float64},D::Vector{Float64},P_r::Vector{Float64},P_i::Vector{Float64},changeling::MMatrix{3, 3, Float64, 9},E_data::Vector{Float64},D_data::Vector{Float64})

  for tii=1:Nt
    update_B(E,B)
    update_uvec(E,D,P_r,P_i,changeling)
    update_E(E,B,P_r,P_i)
    rec = tii/Sampling
    if isinteger(rec) 
      E_data[round(Int,rec)] = E[detection] 
      D_data[round(Int,rec)] = D[detection] 
    end
  end
end


# Run the simulation
t1=time()
#matwrite("D0.mat", Dict("D0" => D0_array)); 
run(E,B,D,P_r,P_i,changeling,E_data,D_data)
t2=time()
println(t2-t1)

matwrite("D.mat", Dict("data" => D_data)); 
matwrite("E.mat", Dict("data" => E_data)); 

matwrite("xE.mat", Dict("data" => E)); 
matwrite("xB.mat", Dict("data" => B)); 
matwrite("xD.mat", Dict("data" => D)); 
matwrite("xP_r.mat", Dict("data" => P_r)); 
matwrite("xP_i.mat", Dict("data" => P_i)); 
