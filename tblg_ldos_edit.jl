using Pkg; Pkg.activate("/Users/mbabar/Desktop/PhD/Analysis/Gerischer/")

loc = "/Users/mbabar/Desktop/PhD/Analysis/Gerischer/ElectrochemicalKinetics.jl/src/"

#using ElectrochemicalKinetics
using MAT
using CSV
using DelimitedFiles
using Printf
using Interpolations
using Trapz

# Compute QC interp
include(loc*"dos.jl")
export DOS
using .DOS: DOSData, get_dos
export DOSData, get_dos
include(loc*"quantum_capacitance.jl")
include(loc*"kinetic_models.jl")
export MarcusHushChidseyDOS

# Testing
#mhcd = MarcusHushChidseyDOS(20, 0.2, "/Users/mbabar/Desktop/PhD/Analysis/Gerischer/ElectrochemicalKinetics.jl/data/DOSes/bernal_graphene.txt")

"""
function calculate_Vdl_interp(dos_f, Vq_min, Vq_max, C_dl)
    Vq_range = range(Vq_min, Vq_max, step=0.001)
    E_min = dos_f.ranges[1][1]
    E_max = dos_f.ranges[1][end]
    CQ = [compute_cq(E_min, E_max, V, dos_f)*1.602*10 for V in Vq_range]
    Vdl_data = []
    Vappl_data = []
    for i = 1:length(Vq_range)
        V_app_i =  Vq_range[i]*((CQ[i]/C_dl) + 1)
        V_dl_i = V_app_i - Vq_range[i]
        push!(Vdl_data, V_dl_i)
        push!(Vappl_data, V_app_i)
    end
    return Vappl_data, Vdl_data
end

function QC_integrand(E, eVq, dos_func; kT=.026)
    dos_func(E) * sech((E-eVq)/(2*kT))^2/(4*kT)
end

function compute_cq(E_min, E_max, eVq, dos_func; kT=.026)
    fcn = E->QC_integrand(E, eVq, dos_func; kT=kT)
    quadgk(fcn, E_min, E_max)[1]
end
"""

Vq_min=-0.45
Vq_max=0.45
C_dl=10.0
#Vappl_data, Vdl_data = calculate_Vdl_interp(mhcd.dos.interp_func, Vq_min, Vq_max, C_dl)
#v_interp = LinearInterpolation(Vappl_data, Vdl_data)

## Debug Integrand limits

function integrand(
    mhcd::MarcusHushChidseyDOS,
    Eo,
    V_dl,
    ox::Bool;
    kT::Real = 0.026,
    V_q = 0.0,
)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -(( E .+ mhcd.λ) .^ 2) ./ (4 * mhcd.λ * kT)
        else
            exp_arg = -(( E .- mhcd.λ) .^ 2) ./ (4 * mhcd.λ * kT)
        end
        exp.(exp_arg)
    end
    fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    E -> mhcd.A .* ((mhcd.dos.interp_func.(E .- V_q)).^ 1) .* marcus_term(E .- (Eo .+ V_q .+ V_dl)) .* fd(E)
end

function compute_k_cq(
    V_app,
    model::MarcusHushChidseyDOS,
    ox::Bool;
    Eo = -0.07, # E_f,red (solvent) - E_f,vac (bilayer) 
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    kT = 0.026,
    E_min = model.dos.E_min,
    E_max = model.dos.E_max,
)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    #Vappl_data, Vdl_data = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    #v_interp = LinearInterpolation(Vappl_data, Vdl_data)

    V_dl = V_dl_interp(V_app)
    #V_dl = v_interp(V_app)
    
    V_q = V_app - V_dl
    if V_q < 0
           E_max = E_max .+ V_q	   
    elseif V_q > 0
    	   E_min = E_min .+ V_q
    end
    k_rate = quadgk(integrand(model, Eo, V_dl, ox; kT = kT, V_q = V_q), E_min, E_max)[1]
    return k_rate, V_q
end

# Modify string
function chop_str(str::String)
         while str[length(str)] == '0'
               str = chop(str)
         end
         if str[length(str)] == '.'
            str = chop(str)
         end
         return str
end

## Calculate rates
# Load matlab ldos data
dos_loc = "/Users/mbabar/Google Drive/My Drive/Analysis/TBLG_Carr/kp_tblg/MATLAB_vers/Data_24kpts/"
theta_list = [0.42, 0.77, 1.15, 1.34, 1.5, 2.39, 2.6, 3, 5];

#theta = theta_list[3];
k_ox = zeros(Float64, size(theta_list));
k_red = zeros(Float64, size(theta_list));
dos_ar = zeros(Float64, size(theta_list)); # ldos area till fermi level

# Interpolate Fermi level
ef_data=[[0.22 -2.7201e-18];[0.43 -0.0017];[0.57 -0.0042];[0.8 -0.0057];[1.1 -0.0078];[1.2 -0.0084];[1.36 -0.0091];[1.55 -0.0098];[1.96 -0.0108];[2.36 -0.0113];[2.67 -0.0116];[2.8 -0.0117];[3.2 -0.0119];[3.73 -0.0121];[4.4 -0.0122];[5 -0.0123]];

ef_func = LinearInterpolation(ef_data[:,1], ef_data[:,2])
kT = 0.026 #eV
lim = 1*kT # limit of dos integration
#print(ef_func(3.4))

## Identify ldos index for AA/AB
region = "AB"
# AA ldos at (1,1) in data matrix
# AB ldos at (12,12) in data matrix
if region=="AA"
   id_x = 1
   id_y = 1
elseif region=="AB"
   id_x = 12
   id_y	= 12
end

Vapp = 0.07 #V
for i in 1:size(theta_list)[1]
    local theta, k, ef, ef_id
    theta = chop_str(string(theta_list[i]));
    file = matopen(dos_loc*"ldos-"*string(theta)*"_90x90.mat")
    data = read(file, "data");
    rscx = read(file, "rscx");
    rscy = read(file, "rscy");

    # Read energy and dos csv
    d = readdlm(dos_loc*"dos-"*string(theta)*".csv", ',', Float64);
    Elist = d[:,1];
    tdos = d[:,2];

    ldos = [Elist data[:,id_x,id_y]]
    ef = ef_func(theta_list[i])
    ef_id = findmin(abs.(ldos[:,1] .- ef))[2] # index closest to fermi level
    lim_id = findmin(abs.(ldos[:,1] .- ef .+ lim))[2] # index closest to integration limit
    
    dos_ar[i] = trapz(ldos[lim_id:ef_id,1], ldos[lim_id:ef_id,2]) # ldos area from ef-kT to ef
    mhcd = MarcusHushChidseyDOS(1.0, 0.82, ldos, Ef=ef)
    k_ox[i], V_q = compute_k_cq(Vapp, mhcd, true; Eo=-0.07,  Vq_min=-0.45, Vq_max=0.45)
    k_red[i], V_q = compute_k_cq(Vapp, mhcd, false; Eo=-0.07,  Vq_min=-0.45, Vq_max=0.45)
    print(i," ",theta," ",k_ox[i]," ",k_red[i]," ", V_q, "\n")
end

#print(dos_ar,"\n","\n")
print(k_ox,"\n")
print(k_red,"\n")
print(k_ox-k_red,"\n")

#print(k)

# Write rate data as .mat
#file = matopen("k_data.mat", "w")
#write(file, "k_data", k_mat)
#write(file, "rscx", rscx)
#write(file, "rscy", rscy)
#close(file)


##
#@time begin
#     mhcd = MarcusHushChidseyDOS(20.0, 0.82, dos)
#     k = compute_k(0.4, mhcd; calc_cq=true, Vq_min=-0.45, Vq_max=0.45)
#end


