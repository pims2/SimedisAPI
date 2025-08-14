################################################################################
# genPatients.jl
#
# Description:
#     Generates simulated patient data for a given scenario by reading threat
#     configurations, initializing the Simedis simulation environment, creating
#     threats, running the simulation, and aggregating victims output files into
#     a single file.
#
# Usage:
#     julia genPatients.jl
#
# Dependencies:
#     CSV, DataFrames, Dates, DelimitedFiles, Logging, OpenStreetMapX, Random,
#     ResumableFunctions, SQLite, SimJulia, TOML, Simedis, Plots, Geodesy, JSON
#
# Author:
#     [Your Name]
# Created:
#     [Date]
################################################################################

# ==============================
# Imports
# ==============================
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

# ==============================
# Global constants & parameters
# ==============================

global victims = []  # Global list to store victim data
consts, paramset = Simedis.constants_from_file("paramsetDemo.toml")  # Load constants and parameters

# ==============================
# Load scenario and threat configs
# ==============================

data = JSON.parsefile("scenario.json")
threat_config_data = JSON.parsefile("threatconfig.json")

# Extract relevant scenario components
threats   = data["threats"]
ccps      = data["ccps"]
hospitals = data["hospitals"]
home      = data["home"]

# Coordinates setup
locations    = [threat["coords"] for threat in threats]
ccp          = [ccp_entry["coords"] for ccp_entry in ccps]
fmp          = [hospital["coords"] for hospital in hospitals]
home_coord   = home
home         = [[50.8467, 4.3525], [50.8467, 4.3525]]  # Default fallback home location

# Convert coordinates to UTM & LLA formats
locationUTM, fmpUTM, ccpUTM, homeUTM, ccp_LLA, fmp_LLA, home_LLA, locations_LLA =
    Simedis.setCoords(locations, ccp, fmp, home, consts.quadrant)

# ==============================
# Map setup
# ==============================
m = Simedis.SetMap("bx.osm")  # Load map file

# ==============================
# Initialize routing & travel times
# ==============================
routes = Vector{Vector{Int}}()
routesH = Vector{Vector{Int}}()

timeCCPHome = []
timeFMPHome = []
timeccpFMP  = []

# Initialize travel times for each threat location
for i in 1:length(locations)
    push!(timeCCPHome, 0.0)
    push!(timeFMPHome, 0.0)
    push!(timeccpFMP, 0.0)
end

# ==============================
# Database setup
# ==============================
baseoutputname = "/home/pims/SimedisAPI/out.sqlite"

# Julia quirk fix for empty type promotion in SQLite operations
Base.promote_type(::Type{Union{}}, ::Type{String}) = String

# Random seed initialization
randomseeds = [abs(rand(Int32)) for i=1:consts.nrruns]
Random.seed!(randomseeds[1])

# Simulation data containers
SimScore_plot = []
SimParams = []
treatment_data = []
global routeG = Vector{Vector{Int}}()

# Input DB connection
inputDB = SQLite.DB(consts.inputDBname)

# ==============================
# Simulation environment setup
# ==============================
local_sim = Simulation()

# Create output database
outputDB = Simedis.createOutputDataBase(consts, baseoutputname)
SQLite.execute(outputDB, "BEGIN TRANSACTION")

# Resource initialization
local_firefighters = Resource(local_sim, consts.nrFF)
local_mediclist    = Simedis.MedicList(local_sim, paramset.Policy, consts)
local_ambulist     = Simedis.AmbuList(local_sim, paramset)
global routesH     = Vector{Vector{Int}}()
patients           = []
local_boatlist     = Simedis.BoatList(local_sim, patients)
local_medevac      = Simedis.Medevac(local_sim)

# Queues for patient transport
local_transportqueue  = []
local_transportqueue1 = []
local_transportqueueT3 = []

# Hospital queue setup
routesH, local_hospitalQueue =
    Simedis.CreateHospitalQueue(inputDB, ccp_LLA, m, routesH, paramset, true)

# ==============================
# Threat creation loop
# ==============================
for i in 1:length(locations)

    # Map textual threat type to numeric code
    subtype = threats[i]["threatType"]
    threat_type_code = if subtype == "artillery strike"
        1
    elseif subtype == "GB release"
        3
    elseif subtype == "drone strike"
        6
    else
        @warn "Unknown threat subtype: $subtype"
        1  # Default fallback
    end

    # Get threat ID from config or use fallback naming
    threat_id = haskey(threats[i], "id") ? string(threats[i]["id"]) : "threat$(i)"

    # Validate threat config presence
    if !haskey(threat_config_data, threat_id)
        error("Missing threat configuration for threat id: $threat_id in threatconfig.json")
    end

    threat_conf = threat_config_data[threat_id]

    # Validate required fields in config
    for key in ["numberOfHits", "scheduledTime"]
        if !haskey(threat_conf, key)
            error("Missing '$key' for threat id $threat_id in threatconfig.json")
        end
    end

    # Extract & convert config values
    number_of_hits = Float64(threat_conf["numberOfHits"])
    scheduled_time = Float64(threat_conf["scheduledTime"])

    # Create threat process in simulation
    @process Simedis.createThreat(
        local_sim,
        threat_type_code,
        scheduled_time,
        locationUTM[i],
        number_of_hits,
        outputDB,
        local_firefighters,
        local_boatlist,
        local_mediclist,
        local_ambulist,
        local_medevac,
        local_transportqueue,
        local_transportqueue1,
        local_transportqueueT3,
        local_hospitalQueue,
        paramset,
        SimScore_plot,
        treatment_data,
        consts,
        inputDB,
        m,
        fmp_LLA[i],
        ccp_LLA[i],
        homeUTM[i],
        ccpUTM[i],
        timeCCPHome[i],
        timeccpFMP[i],
        timeFMPHome[i],
        fmpUTM[i],
        i,
        true
    )
end

# ==============================
# Run the simulation
# ==============================
run(local_sim)

# ==============================
# Victim file aggregation
# ==============================

# Get all victim files matching victims*.txt
victim_files = sort(filter(f -> occursin(r"victims\d+\.txt", f), readdir(".")))

# Identify victims1.txt and victimsN.txt where N > 1
victim1_file     = "victims1.txt"
victim_files_gt1 = filter(f -> parse(Int, match(r"victims(\d+)\.txt", f).captures[1]) > 1, victim_files)

# Aggregate victims into victimsAll.txt
open("victimsAll.txt", "w") do outfile
    # Write victims1.txt content first
    for line in eachline(victim1_file)
        println(outfile, line)
    end

    # Append victims2+.txt content (excluding headers)
    for file in victim_files_gt1
        lines = readlines(file)
        if length(lines) > 1
            for line in lines[2:end]
                println(outfile, line)
            end
        end
    end
end

# Commit simulation results to DB
SQLite.execute(outputDB, "COMMIT")

