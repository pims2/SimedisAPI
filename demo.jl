# Packages
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
#include("ambu.jl")

consts, paramset = Simedis.constants_from_file("paramsetDemo.toml")
#disable_logging(Logging.Warn)

global victims=[]

data = JSON.parsefile("scenario.json")

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
coords=[locations_LLA,fmp_LLA,ccp_LLA,home_LLA]
nSites=length(locations)

inputDB = SQLite.DB(consts.inputDBname)
m = Simedis.SetMap("bx.osm")
routes = Vector{Vector{Int}}()

routes = Vector{Vector{Int}}()
routesH = Vector{Vector{Int}}()

timeCCPHome = []
timeFMPHome = []
timeccpFMP = []
for i in 1:nSites
    push!(timeCCPHome,0.0)
    push!(timeFMPHome,0.0)
    push!(timeccpFMP,0.0)
end

baseoutputname  ="SimBF_$(consts.nrruns)_"*"$(Date(now())).sqlite" # "Output\\000-test"*".sqlite"

Base.promote_type(::Type{Union{}}, ::Type{String}) = String
  
randomseeds = [abs(rand(Int32)) for i=1:consts.nrruns]
Random.seed!(randomseeds[1])

SimScore_plot =[]
#params_data=[]


SimParams=[]
treatment_data = []
global routeG=Vector{Vector{Int}}()
inputDB = SQLite.DB(consts.inputDBname)


outputDB = Simedis.createOutputDataBase(consts,baseoutputname)

SQLite.execute( outputDB, "BEGIN TRANSACTION" )
for nramb in paramset.nrAmbuList
    paramset.nrAmbuFW = nramb  
    paramset.nrAmbuTAC = 0
    for tourniquet in paramset.tourniquetList
        paramset.Tourniquet = tourniquet
    for mascals in paramset.mascalList
            paramset.mascal = mascals
            for Policy1 in paramset.policyList # ["ScoopRun", "StayPlay"]

                if Policy1 == "StayPlay"
                    paramset.Policy = "StayPlay"
                
                elseif Policy1 == "ScoopRun"
                    paramset.Policy = "ScoopRun"
                    
                end

                for PreTriage in paramset.pretriageList # [true, false]
                    paramset.PreTriage = PreTriage
        
                                        for TranspsupervisionLevel in paramset.supervisionList # ["Low", "Medium", "Normal", "High", "Flex"]
                                            paramset.TranspsupervisionLevel = TranspsupervisionLevel

                                            if TranspsupervisionLevel == "Normal"
                                                paramset.supervisionlevels = [[2,1],2,3] # [[2,1],2,3]
                                            elseif TranspsupervisionLevel == "Low"
                                                paramset.supervisionlevels = [3,3,3]
                                            elseif TranspsupervisionLevel == "Medium"
                                                paramset.supervisionlevels = [[2,1],3,3]
                                            elseif TranspsupervisionLevel == "High"
                                                paramset.supervisionlevels = [1,2,3]
                                            elseif TranspsupervisionLevel == "Flex"
                                                paramset.supervisionlevels = [[1,2,3],[2,3],3]
                                            end

                                            for HospitalDistribution in paramset.hospitaldistrList # ["SpreadOut", "CloseFirst"]
                                                paramset.HospitalDistribution = HospitalDistribution

                                                for HospitalCapacity in paramset.hospitalcapalist # ["Low", "Medium", "High"]
                                                    paramset.HospitalCapacity = HospitalCapacity
                                                    
                                                    @time begin
                                                        for k = 1:consts.nrruns

                                                        # if simurun > 2200-1 # custom select simulation run(s)
                                                        # new random random seed for each simulation run
                                                            seed = randomseeds[k] #  j  #abs(rand(Int32))
                                                            
                                                            Random.seed!(seed)
                                                            paramset.seed = seed
                                                        # all simulation code here
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
                                                            if paramset.simurun==1
                                                                global routeG = routesH
                                                            end
                                                            paramset.hospitalOrder = [hospital.ID for hospital in local_hospitalQueue]

                                                            if !consts.keepAllInfo
                                                                SQLite.execute( outputDB, "DELETE FROM VictimLog")
                                                                SQLite.execute( outputDB, "DELETE FROM ResourceLog")
                                                            end
                                                            if consts.creator 
                                                                @process Simedis.createThreat(local_sim,5,0.0,locationUTM[1],1.0,outputDB,local_firefighters,local_boatlist, local_mediclist, local_ambulist,local_medevac, local_transportqueue, local_transportqueue1,local_transportqueueT3, local_hospitalQueue, paramset, SimScore_plot,treatment_data,consts,inputDB,m,fmp_LLA[1],ccp_LLA[1],homeUTM[1],ccpUTM[1],timeCCPHome[1],timeccpFMP[1],timeFMPHome[1],fmpUTM[1],0,true) 
                                                                    
                                                            #    @process Simedis.createThreat(local_sim,6,480.0,locationUTM[2],1.0,outputDB,local_firefighters,local_boatlist, local_mediclist, local_ambulist,local_medevac, local_transportqueue,local_transportqueue1, local_transportqueueT3, local_hospitalQueue, paramset, SimScore_plot,treatment_data,consts,inputDB,m,fmp_LLA[2],ccp_LLA[2],homeUTM[2],ccpUTM[2],timeCCPHome[2],timeccpFMP[2],timeFMPHome[2],fmpUTM[2],1,true)
                                                            
                                                                
                                                            else
                                                                @process Simedis.readVictimsMulti(local_sim, local_firefighters, local_mediclist, local_ambulist,local_medevac, local_transportqueue, local_transportqueueT3, local_hospitalQueue, paramset,SimScore_plot,treatment_data,consts,outputDB,fmpUTM,homeUTM,ccpUTM,fmp_LLA,ccp_LLA,timeCCPHome,timeccpFMP,timeFMPHome,inputDB,m,5,locations,true)
                                                            end
                                                            
                                                            @process Simedis.scheduleAll(local_sim, local_firefighters, local_mediclist, local_boatlist,local_ambulist,local_medevac, local_transportqueue, local_transportqueueT3, local_hospitalQueue, paramset, inputDB,consts,treatment_data,outputDB,ccp,m,false,false,homeUTM)
                                                            
                                                            run(local_sim)

                                                        # custom select simulation run(s)
                                                    
                                                        Simedis.AddOutputOverviewDB(paramset, outputDB, local_mediclist, local_ambulist,consts)
                                                        Simedis.reportLocations(inputDB,outputDB,m,paramset,1000,routeG,"bxl",coords,consts,0.0,0.0,true,true) # plot the html file
                                                        
                                                        
                                                    
                                                        if mod(paramset.simurun,100) == 0
                                                            SQLite.execute( outputDB, "COMMIT" )
                                                            SQLite.execute( outputDB, "BEGIN TRANSACTION" )
                                                        end
                                                
                                                            paramset.simurun += 1
                                                            
                                                        #Output Files
                                                        if consts.writeSimScore
                                                            writedlm("simedisscoredata.txt",SimScore_plot)
                                                        end
                                                        #writedlm("treatment_data.txt",treatment_data)
                                                        # writedlm("params.txt",params_data)
                                                    
                                                        # writedlm("triages.txt")
                                                        # if consts.creator
                                                        #     writedlm("victims.txt",victims_list)
                                                        # end
                                                        GC.gc()
                                                        # nrruns
                                                        end
                                                    end
                                                end
                                                #hospitalCapacity
                                            end
                                            # HospitalDistribution
                                        end
                                        # Supervision
                                    
                end
                # PreTriage
        end
        # Policy
    end
    #Mascal
    end
    #tourniquet
end

SQLite.execute( outputDB, "COMMIT" )
println("the end....")