using CSV
using DataFrames
using Dates
using DelimitedFiles
using Logging
using OpenStreetMapX
using Random
using ResumableFunctions
using SQLite
using SimJulia
using TOML
using Simedis
using Plots
using Geodesy
using JSON


global victims=[]
consts, paramset = Simedis.constants_from_file("paramsetDemo.toml")


data =JSON.parsefile("scenario.json")

threats = data["threats"]
ccps = data["ccps"]
hospitals = data["hospitals"]
home = data["home"]

locations = [threat["coords"] for threat in threats]
ccp = [ccp_entry["coords"] for ccp_entry in ccps]
fmp = [hospital["coords"] for hospital in hospitals]
home_coord = home
home=[[50.8467, 4.3525],[50.8467, 4.3525]]
locationUTM,fmpUTM,ccpUTM,homeUTM,ccp_LLA,fmp_LLA,home_LLA,locations_LLA= Simedis.setCoords(locations,ccp,fmp,home,consts.quadrant)
println(locationUTM,fmpUTM,ccpUTM,homeUTM,ccp_LLA,fmp_LLA,home_LLA,locations_LLA)
m=Simedis.SetMap("demo.osm") #peu importe


routes = Vector{Vector{Int}}()
routesH = Vector{Vector{Int}}()

timeCCPHome = []
timeFMPHome = []
timeccpFMP = []
for i in 1:length(locations)
    push!(timeCCPHome,0.0)
    push!(timeFMPHome,0.0)
    push!(timeccpFMP,0.0)
end
# coords=[locations_LLA,fmp_LLA,ccp_LLA,home_LLA]
baseoutputname  ="/home/pims/simedis-ui/runSimedis/out.sqlite" # "Output\\000-test"*".sqlite"

Base.promote_type(::Type{Union{}}, ::Type{String}) = String
  
randomseeds = [abs(rand(Int32)) for i=1:consts.nrruns]
Random.seed!(randomseeds[1])

SimScore_plot =[]
#params_data=[]


SimParams=[]
treatment_data = []
global routeG=Vector{Vector{Int}}()
inputDB = SQLite.DB( consts.inputDBname )

outputDB = Simedis.createOutputDataBase(consts,baseoutputname)

SQLite.execute( outputDB, "BEGIN TRANSACTION" )
local_sim = Simulation()
local_firefighters = Resource(local_sim, consts.nrFF)
local_mediclist = Simedis.MedicList(local_sim, paramset.Policy,consts)
local_ambulist = Simedis.AmbuList(local_sim,paramset)
global routesH=Vector{Vector{Int}}()
patients=[]
local_boatlist = Simedis.BoatList(local_sim,patients)
local_medevac = Simedis.Medevac(local_sim)
local_transportqueue = []
local_transportqueue1 = []
local_transportqueueT3 = []
routesH,local_hospitalQueue = Simedis.CreateHospitalQueue(inputDB, ccp_LLA,m,routesH,paramset,true)
for i in 1:length(locations)
    subtype = threats[i]["threatType"]
    threat_type_code = 0
    if subtype == "artillery strike"
        threat_type_code = 1
    elseif subtype == "GB release"  # GB release
        threat_type_code = 3
    elseif subtype == "drone strike"
        threat_type_code = 6
    else
        @warn "Unknown threat subtype: $subtype"
        threat_type_code = 1  # default
    end 
    if threat_type_code == 1
        @process Simedis.createThreat(local_sim,threat_type_code,0.0,locationUTM[i],3.0,outputDB,local_firefighters,local_boatlist, local_mediclist, local_ambulist,local_medevac, local_transportqueue, local_transportqueue1,local_transportqueueT3, local_hospitalQueue, paramset, SimScore_plot,treatment_data,consts,inputDB,m,fmp_LLA[i],ccp_LLA[i],homeUTM[i],ccpUTM[i],timeCCPHome[i],timeccpFMP[i],timeFMPHome[i],fmpUTM[i],i,true) 
    else
         @process Simedis.createThreat(local_sim,threat_type_code,0.0,locationUTM[i],1.0,outputDB,local_firefighters,local_boatlist, local_mediclist, local_ambulist,local_medevac, local_transportqueue, local_transportqueue1,local_transportqueueT3, local_hospitalQueue, paramset, SimScore_plot,treatment_data,consts,inputDB,m,fmp_LLA[i],ccp_LLA[i],homeUTM[i],ccpUTM[i],timeCCPHome[i],timeccpFMP[i],timeFMPHome[i],fmpUTM[i],i,true) 
    end
end
run(local_sim)
SQLite.execute( outputDB, "COMMIT" )
