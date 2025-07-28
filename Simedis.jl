module Simedis

using CSV
using Configurations
using DataFrames
using Dates
using DelimitedFiles
using Distributions
using Geodesy
using Logging
using OpenStreetMapX
using Profile
using Random
using ResumableFunctions
using SQLite
using SimJulia
using SpecialFunctions
using TOML
using TimerOutputs
using Graphs
using Parameters
using PyCall
using Conda
using JSON


@option "Parameterset" mutable struct Parameterset
    Policy :: String
    PreTriage :: Bool
    PreTriageLocation :: String
    TriageLevel :: String
    HospitalDistribution :: String
    HospitalCapacity :: String
    TranspsupervisionLevel :: String
    hospitalOrder :: Array
    Tourniquet :: Bool
    supervisionlevels :: Array
    mascal :: Bool
    buildFMP :: Bool
    nrAmbuFW :: Int
    nrAmbuTAC :: Int
    policyList :: Array
    pretriageList :: Array
    supervisionList :: Array
    hospitaldistrList :: Array
    hospitalcapalist ::Array
    tourniquetList :: Array
    mascalList :: Array
    nrAmbuList :: Array

    seed :: Int64 
    simurun :: Int64 = 1
    TriageTimeUrgent :: Real = 2.0
    
end
@option "Constants" mutable struct Constants
    nrruns :: Int
    numberofVictims :: Int
    creator :: Bool
    writeSimScore :: Bool
    tabularIP :: Bool
    logg :: Bool
    walkSpeed :: Float64
    quadrant :: Int
    offroad :: Bool
    stochasTimes :: Bool
    keepAllInfo :: Bool
    inputDBname :: String
    nrBoats :: Int
    BoatCapacity :: Int
    nrAmbuPreliminary :: Int
    nrFF :: Int
    minNurseFMP :: Int
    minDoctorFMP :: Int
    nrTriageMMT :: Int
    setupFMPtime :: Float64
    LoadFMP :: Float64
    UnloadFMP :: Float64
    DriveTimeFMP :: Float64 #A update par un calcul avec la carte!
    UnloadHosp :: Float64
    HospDropOff :: Float64
    TriageTimeT3 :: Float64
    PreTriageTime :: Float64
    TriageTimeUrgentDefault :: Float64
    sigmaLN :: Float64
    sigmaN :: Int
    TreatTimeDist :: String
    TreatTimeBound :: Float64
    TravTimeBound :: Float64
    TravTimeDist :: String
    RescueTimeBound :: Float64
    RescueTimeDist :: String
    OtherTimeDist :: String
    OtherTimeBound :: Float64
    disrobetimes :: Array
    passtimes :: Array
    decontimes ::Array
    rerobetimes ::Array
    amst_arrivaltime :: Float64
    ambushRate :: Float64
         
end

mutable struct Event
    name :: String #definition of the event
    aoe :: Float64 #area of effect in squared meters i.e. blast radius
    scheduledTime :: Float64 #when is the event programmed?
    magnitude :: Int64 #severity of the event, directly impacts casualty numbers
    location :: UTM #event position (ex. blast center)
    cbrn :: Int64 #if there is a CBRN threat define it here
    threatLength :: Float64 #minutes duration of the threat
    sustained :: Bool #is the event sustained, i.e. gunfire where every X minutes casualties are expected
    casualtyRate :: Int64 #estimated #of new casualties for a sustained threat
    casualtyEstimate :: Int64 #estimated total casualtyNumber to generate

    function Event(name:: String, aoe::Float64, scheduledTime::Float64, magnitude ::Int64,location :: UTM, cbrn::Int64, threatLength ::Float64,sustained::Bool, casualtyRate::Int64,consts::Constants)

        event = new()
        event.name = name
        event.aoe = aoe
        event.scheduledTime = scheduledTime
        event.magnitude = magnitude
        event.location = location
        event.cbrn = cbrn #1 for Chemical, 2 for Bacteriological, 3 radiological, 4 Nuclear, 5 Explosion, 0 for none
        event.threatLength = threatLength
        event.sustained = sustained
        event.casualtyRate= casualtyRate
        event.casualtyEstimate= casualtyRate*threatLength

        return event
    end
end

mutable struct Amst
    atropine ::Int64
    oxygen ::Int64
    
    atropineDoses ::Int64
    capacityOxygen ::Int64
    atropineInjectionReady :: Bool
    
    function Amst(env ::Environment, atropineDoses ::Int64,  capacityOxygen ::Int64,atropineInjectionReady::Bool)
        ams = new()
        ams.atropine = atropineDoses
        ams.oxygen = capacityOxygen
        ams.atropineInjectionReady=atropineInjectionReady
        return ams
    end
end
mutable struct DeconTrack
    disrobe ::Container
    decon ::Container
    rerobe ::Container
    capacitydisrobe ::Int64
    capacitydecon ::Int64
    capacityrerobe ::Int64

    function DeconTrack(env ::Environment, capacitydisrobe ::Int64,  capacitydecon ::Int64, capacityrerobe ::Int64)
        track = new()
        track.disrobe = Resource(env,capacitydisrobe)
        track.decon = Resource(env,capacitydecon)
        track.rerobe = Resource(env,capacityrerobe)
        return track
    end
end
mutable struct DeconUnit
    disrobeAvailableIncap ::Container
    decontracksIncap ::Array{DeconTrack,1}
    disrobeAvailableMobile ::Container
    decontracksMobile ::Array{DeconTrack,1}
    nrtracksIncap ::Int64
    capacitydisrobeIncap :: Int64
    capacitydeconIncap ::Int64
    capacityrerobeIncap ::Int64
    nrtracksMobile ::Int64
    capacitydisrobeMobile ::Int64
    capacitydeconMobile ::Int64
    capacityrerobeMobile ::Int64

    function DeconUnit(env ::Environment, nrtracksIncap ::Int64, capacitydisrobeIncap ::Int64,  capacitydeconIncap ::Int64, capacityrerobeIncap ::Int64, nrtracksMobile ::Int64, capacitydisrobeMobile ::Int64,  capacitydeconMobile ::Int64, capacityrerobeMobile ::Int64)
        unit = new()
        unit.capacitydeconMobile = capacitydeconMobile
        unit.capacitydeconIncap = capacitydisrobeIncap
        unit.capacitydisrobeIncap = capacitydisrobeIncap
        unit.capacitydisrobeMobile = capacitydisrobeMobile
        unit.capacityrerobeIncap = capacityrerobeIncap
        unit.capacityrerobeMobile = capacityrerobeMobile

        if nrtracksIncap < 1
            unit.disrobeAvailableIncap = Resource(env,0)
            unit.decontracksIncap = []
        else
            unit.disrobeAvailableIncap = Resource(env,1)
            for i in 1:nrtracksIncap
                if i == 1
                    unit.decontracksIncap = [DeconTrack(env, capacitydisrobeIncap,  capacitydeconIncap, capacityrerobeIncap)]
                else
                    push!(unit.decontracksIncap, DeconTrack(env, capacitydisrobeIncap,  capacitydeconIncap, capacityrerobeIncap))
                end
            end
        end
        if nrtracksMobile < 1
            unit.disrobeAvailableMobile = Resource(env,0)
            unit.decontracksMobile = []
        else
            unit.disrobeAvailableMobile = Resource(env,1)
            for i in 1:nrtracksMobile
                if i == 1
                    unit.decontracksMobile = [DeconTrack(env, capacitydisrobeMobile,  capacitydeconMobile, capacityrerobeMobile)]
                else
                    push!(unit.decontracksMobile, DeconTrack(env, capacitydisrobeMobile,  capacitydeconMobile, capacityrerobeMobile))
                end
            end
        end

        return unit
    end
end

mutable struct Injury
    tag :: Int64 #trauma tag from list of injuries
    injuryLocation :: String #location of the injury on the body
    name :: String #definition of the injury
    severity :: String #minor/moderate/severe
    severityScore :: Int64 #injury Abbreviated Injury Score
    role::Array
    treatmentLocation :: String #place where the treatment can be administered
    treatmentTime :: Float64
    hospitalNeed :: String #hospital type needed depending on injury
    equipmentMitigation :: Float64
    
    
    function Injury(tag ::Int64,injuryLocation ::String,name ::String,severity ::String,severityScore :: Int64,role :: Array,treatmentLocation::String,treatmentTime::Float64,hospitalNeed::String,equipmentMitigation::Float64)

        injury = new()
        injury.tag = tag
        injury.injuryLocation = injuryLocation
        injury.name = name
        injury.severity=severity
        injury.severityScore=severityScore
        injury.role=role
        injury.treatmentLocation = treatmentLocation
        injury.treatmentTime = treatmentTime
        injury.hospitalNeed = hospitalNeed
        injury.equipmentMitigation = equipmentMitigation
       
        return injury
    end

end
mutable struct Patient

    env :: Environment
    vicid ::Int64 #tag to name a victim, unique
    age :: Int64
    is_treated :: Int64
    mobile :: Int64
    isLethallyInjured :: Int64
    triage :: Int64
    x0 :: Float64
    y0 :: Float64
    iss :: Int64
    timeOfDeath :: Float64
    scheduledTime :: Float64
    SimParams :: Array
    maxSimScore :: Float64
    isBleeding :: Int64
    timetrigtransit ::Process
    effecttimeproc ::Process
    transptreatmentproc ::Process
    transport ::Container
    treatmenthasstarted ::Bool
    TriageIncorrectProb :: Tuple
    triagecorrectness ::String
    tourniquet :: Container
    combopen :: Container
    victimInjury :: Injury
    victimInjury2 :: Injury
    victimInjury3 :: Injury
    hospitalneed :: String
    hospitalID :: Int64
    bloodvolume :: Float64
    incapacitated :: Int64
    treatmentEndTime :: Float64
    isdead :: Int64 
    hasCrossed :: Bool
    SimScore_plot :: Array
    distFromBlast :: Float64
    treatmentParams :: Array
    disasterSite :: Int64
    facility :: String
    dose :: Float64
    ip :: Int64
    contaminated :: Bool
   
        
   
    function Patient(env :: Environment, vicid :: Int64, age :: Int64, is_treated::Int64,   mobile ::Int64,  isLethallyInjured:: Int64, triage:: Int64 ,x0 :: Float64,y0 :: Float64,iss::Int64,timeOfDeath::Float64,scheduledTime::Float64,SimParams :: Array,maxSimScore::Float64,isBleeding::Int64, TriageIncorrectProb :: Tuple,victimInjury1::Injury,victimInjury2::Injury,victimInjury3::Injury,hospitalneed::String,consts::Constants,SimScore_plot,distFromBlast,treatmentParams,disasterSite::Int64,facility::String,ip::Int64,contaminated::Bool=false,dose::Float64=0.0,hospitalID::Int64=0, bloodvolume::Float64=5.0,incapacitated::Int64=0,treatmentEndTime::Float64=0.0, isdead::Int64 = 0,hasCrossed::Bool =false)
        patient = new() # create new patient, with all properties empty
        patient.vicid = vicid
        patient.age = age
        patient.is_treated = is_treated
        patient.mobile = mobile #default patient is able to walk
        patient.isLethallyInjured = isLethallyInjured
        patient.triage = triage
        patient.x0 = x0
        patient.y0 = y0
        patient.iss = iss
        patient.timeOfDeath = timeOfDeath
        patient.scheduledTime = scheduledTime #time of inject in minutes offset from start
        patient.maxSimScore = maxSimScore
        patient.isBleeding = isBleeding
        patient.transport = Resource(env,1)
        patient.treatmenthasstarted = false
        patient.TriageIncorrectProb = (0,0)
        patient.triagecorrectness = DetermineTriageCorrectness(TriageIncorrectProb)
        patient.tourniquet = Resource(env,1)  #each soldier has a tourniquet
        patient.combopen = Resource(env,1) #each soldier has a combopen
        patient.victimInjury = victimInjury1
        patient.victimInjury2 = victimInjury2
        patient.victimInjury3 = victimInjury3
        patient.distFromBlast = distFromBlast
        patient.treatmentParams = treatmentParams
        patient.dose=dose
        patient.distFromBlast = distFromBlast
        patient.treatmentParams = treatmentParams
        patient.dose=dose
        patient.hospitalneed = hospitalneed
        patient.hospitalID = 0
        patient.bloodvolume = bloodvolume
        patient.facility = facility
        patient.ip=ip
        patient.contaminated=contaminated
        patient.incapacitated = incapacitated
        patient.treatmentEndTime = treatmentEndTime
        patient.isdead = isdead
        patient.hasCrossed = hasCrossed
        patient.disasterSite = disasterSite #from which incident
        if consts.creator ===true
             SimParams = clinconds(env, patient,SimScore_plot,consts) #defines the SimedisScore evolution of the patient
             patient.SimParams = SimParams
        else
            patient.SimParams = SimParams
            if consts.writeSimScore
                clinconds(env, patient,SimScore_plot,consts)
            end
            
        end

      
        patient.hasCrossed = hasCrossed
        patient.disasterSite = disasterSite #from which incident
        if consts.creator ===true
             SimParams = clinconds(env, patient,SimScore_plot,consts) #defines the SimedisScore evolution of the patient
             patient.SimParams = SimParams
        else
            patient.SimParams = SimParams
            if consts.writeSimScore
                clinconds(env, patient,SimScore_plot,consts)
            end
            
        end

      
        
      
        return patient
    end
end
mutable struct MedicList
        onSiteDoctor :: Container
        onSiteNurse :: Container
        ccpDoctor :: Container
        ccpNurse :: Container
        fmpMMT :: Container
        fmpDoctor :: Container
        fmpNurse :: Container
        nrTriageMMTArrived ::Int
        nrDoctorTransport :: Int
        nrNurseTransport :: Int
        minDoctorFMP ::Int
        minNurseFMP ::Int
        dirMed :: Container
        #assDirMed :: Container
        BuildingFMP :: Bool
        TreatmentAtFMP :: Bool
        

    function MedicList(env :: Environment, policy :: String,consts::Constants)
        list = new()
        list.onSiteDoctor = Resource(env,0)
        list.onSiteNurse = Resource(env,0)
        list.ccpDoctor = Resource(env,0)
        list.ccpNurse = Resource(env,0)
        list.fmpMMT = Resource(env,0)
        list.fmpDoctor = Resource(env,2)
        list.fmpNurse = Resource(env,4)
        list.nrTriageMMTArrived = 0
        list.nrDoctorTransport = 0
        list.nrNurseTransport = 0
        list.minDoctorFMP = consts.minDoctorFMP
        list.minNurseFMP = consts.minNurseFMP
        list.dirMed = Resource(env,1)
        #list.assDirMed = Resource(env,0)

        list.BuildingFMP = false
        if policy == "ScoopRun"
            list.TreatmentAtFMP = false
        elseif policy == "StayPlay"
            list.TreatmentAtFMP = false
        end
        return list
    end
end
mutable struct AmbuList #total pool of ambulances
   
    ForwardMedevac ::Container
    TacMedevac ::Container
    transportchecktrigger ::Container
    transportcheckoccupied ::Container
    
    allready ::Bool #signal that transport may be finished (this is the mission_done equivalent)


    function AmbuList(env :: Environment,paramset::Parameterset)
        list = new()
        list.ForwardMedevac = Resource(env,paramset.nrAmbuFW)
        list.TacMedevac = Resource(env,paramset.nrAmbuTAC)
        list.transportchecktrigger = Resource(env, 1)
        list.transportcheckoccupied = Resource(env, 1)
        list.allready = false
        
        request(list.transportchecktrigger)
        # list.transportcheckoccupied should be available

        return list
    end
end
mutable struct BoatList #total pool of ambulances
   
    ForwardMedevac ::Container
    TacMedevac ::Container
    transportchecktrigger ::Container
    transportcheckoccupied ::Container
    patients :: Array
    
    
    allready ::Bool #signal that transport may be finished (this is the mission_done equivalent)


    function BoatList(env :: Environment,patients::Array)
        list = new()
        list.ForwardMedevac = Resource(env,0)
        list.TacMedevac = Resource(env,0)
        list.transportchecktrigger = Resource(env, 1)
        list.transportcheckoccupied = Resource(env, 1)
        list.patients = patients
        list.allready = false
        
        request(list.transportchecktrigger)
        # list.transportcheckoccupied should be available

        return list
    end
end
mutable struct TruckList #total pool of ambulances
   
    ForwardMedevac ::Container
    TacMedevac ::Container
    transportchecktrigger ::Container
    transportcheckoccupied ::Container
    
    allready ::Bool #signal that transport may be finished (this is the mission_done equivalent)

    function TruckList(env :: Environment)
        list = new()
        list.ForwardMedevac = Resource(env,0)
        list.TacMedevac = Resource(env,0)
        list.transportchecktrigger = Resource(env, 1)
        list.transportcheckoccupied = Resource(env, 1)
        list.allready = false
        
        request(list.transportchecktrigger)
        # list.transportcheckoccupied should be available

        return list
    end
end
mutable struct Medevac
    preliminary::Container
    primary::Container

    transportchecktrigger ::Container
    transportcheckoccupied ::Container
    allready ::Bool

    function Medevac(env::Environment)
        list = new()
        list.preliminary = Resource(env,0)
        list.primary = Resource(env,0)
        list.transportchecktrigger = Resource(env, 1)
        list.transportcheckoccupied = Resource(env, 1)
        list.allready = false
        list.preliminary.capacity = 2
        list.primary.capacity = 2

        request(list.transportchecktrigger)

        return list
    end
end
mutable struct Hospital
    ID :: Int64
    TravelTime :: Float64
    Surgemode :: Bool

    CapT1 :: Int64
    CapT2 :: Int64
    CapT3 :: Int64
    SurgeT1 :: Int64
    SurgeT2 :: Int64
    SurgeT3 :: Int64

    CapT1Ped :: Int64
    CapT2Ped :: Int64
    CapT3Ped :: Int64
    SurgeT1Ped :: Int64
    SurgeT2Ped :: Int64
    SurgeT3Ped :: Int64

    AdmittedT1 :: Int64
    AdmittedT2 :: Int64
    AdmittedT3 :: Int64
    AdmittedT1Ped :: Int64
    AdmittedT2Ped :: Int64
    AdmittedT3Ped :: Int64
    LAT:: Float64
    LON:: Float64
  

    function Hospital(dfHospital :: DataFrameRow{DataFrame,DataFrames.Index}, HospitalCapacity :: String)
        hospital = new()
        hospital.ID = dfHospital.ID
        hospital.TravelTime = dfHospital.DriveTime
        hospital.LAT = dfHospital.LAT
        hospital.LON = dfHospital.LON
       
        hospital.Surgemode = false
        if HospitalCapacity == "Low"

            hospital.CapT1 = dfHospital.CapT1_A_Low
            hospital.CapT2 = dfHospital.CapT2_A_Low
            hospital.CapT3 = dfHospital.CapT3_A_Low

            hospital.CapT1Ped = dfHospital.CapT1_P_Low
            hospital.CapT2Ped = dfHospital.CapT2_P_Low
            hospital.CapT3Ped = dfHospital.CapT3_P_Low

        elseif HospitalCapacity == "Medium"

            hospital.CapT1 = dfHospital.CapT1_A_Medium
            hospital.CapT2 = dfHospital.CapT2_A_Medium
            hospital.CapT3 = dfHospital.CapT3_A_Medium

            hospital.CapT1Ped = dfHospital.CapT1_P_Medium
            hospital.CapT2Ped = dfHospital.CapT2_P_Medium
            hospital.CapT3Ped = dfHospital.CapT3_P_Medium

        elseif HospitalCapacity == "High"

            hospital.CapT1 = dfHospital.CapT1_A_High
            hospital.CapT2 = dfHospital.CapT2_A_High
            hospital.CapT3 = dfHospital.CapT3_A_High

            hospital.CapT1Ped = dfHospital.CapT1_P_High
            hospital.CapT2Ped = dfHospital.CapT2_P_High
            hospital.CapT3Ped = dfHospital.CapT3_P_High

        end

        hospital.SurgeT1 = dfHospital.SurgeCapT1_A
        hospital.SurgeT2 = dfHospital.SurgeCapT2_A
        hospital.SurgeT3 = dfHospital.SurgeCapT3_A

        hospital.SurgeT1Ped = dfHospital.SurgeCapT1_P
        hospital.SurgeT2Ped = dfHospital.SurgeCapT2_P
        hospital.SurgeT3Ped = dfHospital.SurgeCapT3_P

        hospital.AdmittedT1 = 0
        hospital.AdmittedT2 = 0
        hospital.AdmittedT3 = 0
        hospital.AdmittedT1Ped = 0
        hospital.AdmittedT2Ped = 0
        hospital.AdmittedT3Ped = 0
        

        return hospital
    end
end

function constants_from_file(path::AbstractString)
    parsed_toml=TOML.parsefile(path)
    paramset = from_dict(Parameterset, parsed_toml["Parameterset"])
    consts = from_dict(Constants, parsed_toml["Constants"])

    return consts, paramset
end

"""
setCoords(locations::Array,ccp::Array,fmp::Array,home::Array,quadrant::Int64)

Set map coordinates


"""

function setCoords(locations::Array,ccp::Array,fmp::Array,home::Array,quadrant::Int64)

        locations_LLA = []
        ccp_LLA = []
        fmp_LLA = []
        home_LLA = []
        locationUTM = []
        fmpUTM = []
        ccpUTM =[]
        homeUTM = []
       
        
        utm_location = UTMfromLLA(quadrant,true,wgs84)

             
        for i in 1:length(locations)
            location1_lla = Geodesy.LLA(locations[i][1],locations[i][2])
            location1_utm = utm_location(location1_lla)
            
            ccp1_LLA =  Geodesy.LLA(ccp[i][1],ccp[i][2])
            ccp1_utm =utm_location(ccp1_LLA)
            fmp1_LLA = Geodesy.LLA(fmp[i][1],fmp[i][2])
            fmp1_utm = utm_location(fmp1_LLA)
            home1_LLA = Geodesy.LLA(home[i][1],home[i][2])
            home1_utm = utm_location(home1_LLA)
            push!(locations_LLA,location1_lla)
            push!(locationUTM,location1_utm)
            push!(ccp_LLA,ccp1_LLA)
            push!(ccpUTM,ccp1_utm)
            push!(fmp_LLA,fmp1_LLA)
            push!(fmpUTM,fmp1_utm)
            push!(home_LLA,home1_LLA)
            push!(homeUTM,home1_utm)
        end
        
        
    return locationUTM,fmpUTM,ccpUTM,homeUTM,ccp_LLA,fmp_LLA,home_LLA,locations_LLA

end

function SetMap(file::String)
    map=get_map_data(file,use_cache=false, trim_to_connected_graph=true);
    return map

end


function computeBloodVolume(patient::Patient,t::Float64,t_TQ::Float64,alpha::Float64,k::Float64,beta::Float64, controlled::Bool)
    bv0 = 5.0 #%bv in liters
    if controlled
        bv=bv0*exp(-beta*(t-t_TQ))
    else #uncontrolled
        bv=bv0*exp(-(alpha+k*t)*t)
    end

        return bv
end

@resumable function signalDisrobeAvailable(env ::Environment, disrobeAvailable ::Container)
    if disrobeAvailable.level == disrobeAvailable.capacity
        
        @yield release(disrobeAvailable)
    end
end
function findavailabletrack(unit ::DeconUnit, vicmobile ::Bool)
    if vicmobile
        for i in 1:length(unit.decontracksMobile)
            if unit.decontracksMobile[i].disrobe.level < unit.decontracksMobile[i].disrobe.capacity
                return i
            end
        end
    else
        for i in 1:length(unit.decontracksIncap)
            if unit.decontracksIncap[i].disrobe.level < unit.decontracksIncap[i].disrobe.capacity
                return i
            end
        end
    end
    println("No track available")
end
@resumable function wetdecontaminationprocess(env ::Environment, patient ::Patient, unit ::DeconUnit, mobilerequest ::Bool,params,consts,outputDB)

        

        if checkMobility(env,patient)
            disrobesignal = unit.disrobeAvailableMobile
        else
            disrobesignal = unit.disrobeAvailableIncap
        end

        tracknr = findavailabletrack(unit, checkMobility(env,patient))
   

        if mobilerequest
            deconTrack = unit.decontracksMobile[tracknr]
        else
            deconTrack = unit.decontracksIncap[tracknr]
        end

        # Decontamination procedure
        # Disrobe
        @yield request(deconTrack.disrobe)
        # see if other track may be available, if so, signaldisrobe avialable
        try findavailabletrack(unit, mobilerequest)
            @process signalDisrobeAvailable(env, disrobesignal)
        catch
            println("No other tracks")
        end
        
      
        logData(env,patient,8,params,consts,outputDB)
        
        @yield timeout(env, timerandomised(consts.disrobetimes[patient.mobile+1],"Normal",consts))

        logData(env,patient,9,params,consts,outputDB)
        
        # Request decon and release disrobe afterwards
        @yield request(deconTrack.decon)
        @yield timeout(env, timerandomised(consts.passtimes[patient.mobile+1],"Normal",consts))
        @yield release(deconTrack.disrobe)

        @process signalDisrobeAvailable(env, disrobesignal)
   
        # Start WET decontamination
    
        unit.capacitydeconMobile = unit.capacitydeconMobile -1 
      
        logData(env,patient,10,params,consts,outputDB)
        

        @yield timeout(env, timerandomised(consts.decontimes[patient.mobile+1],"Normal",consts))

        logData(env,patient,11,params,consts,outputDB)
                
        # Rerobe
        @yield request(deconTrack.rerobe)
        @yield timeout(env, timerandomised(consts.passtimes[patient.mobile+1],"Normal",consts))
        @yield release(deconTrack.decon)
    
        logData(env,patient,12,params,consts,outputDB)
                    
        @yield timeout(env, timerandomised(consts.rerobetimes[patient.mobile+1],"Normal",consts))
      
        logData(env,patient,13,params,consts,outputDB)
                    
        @yield release(deconTrack.rerobe)

        #patient.decontaminated = true
        #patient.SimParams[5] = 7
   
        
end
@resumable function decontamination(env ::Environment, patient ::Patient, unit ::DeconUnit,consts,params,outputDB)
    # Function that does the decontamination of a victim and possibly the stabilisation on site
    if isdead(patient)
        return
    end
    if checkMobility(env,patient)
        mobilerequest = true
        disrobesignal = unit.disrobeAvailableMobile
    else
        mobilerequest = false
        disrobesignal = unit.disrobeAvailableIncap
    end

    @yield requestdecon = request(disrobesignal; priority = patient.triage)

    wetdeconprocess = @process wetdecontaminationprocess(env, patient, unit, mobilerequest,params,consts,outputDB)
    @yield wetdeconprocess


end
function milInjurySetter(injuryNumber::Int64,victimInjury::Injury)
    #reads an injury number and sets the injury to the victim
    if injuryNumber == 0
        victimInjury.name = "none"
        victimInjury.severity = "none"
        victimInjury.severityScore = 0
        victimInjury.role = ["none"]
        victimInjury.treatmentLocation = "home"
        victimInjury.injuryLocation = "none"
        victimInjury.hospitalNeed = "H"
        victimInjury.equipmentMitigation = 1.0
    end
    if injuryNumber == 1
        victimInjury.name = "shrapnel penetration arm"
        victimInjury.severity = "moderate"
        victimInjury.severityScore = 3
        victimInjury.role = ["role1","role2","role3","role3"]
        victimInjury.treatmentLocation = "ED"
        victimInjury.injuryLocation = "upper extremity"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.8
        
    end
    if injuryNumber == 2
        victimInjury.name = "tympanic membrane rupture"
        victimInjury.severity = "minor"
        victimInjury.severityScore = 1
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "ANY"
        victimInjury.injuryLocation = "head/neck"
        victimInjury.hospitalNeed = "OPC"
        victimInjury.equipmentMitigation = 0.95 #headset+helmet
    end
    
    if injuryNumber == 3
        victimInjury.name = "concussion"
        victimInjury.severity = "minor"
        victimInjury.severityScore = 1
        victimInjury.role = ["analgesics"]
        victimInjury.treatmentLocation = "ANY"
        victimInjury.injuryLocation = "head/neck"
        victimInjury.hospitalNeed = "OPC"
        victimInjury.equipmentMitigation = 0.8 #helmet
    end
    if injuryNumber == 4
        victimInjury.name = "blast lung injury"
        victimInjury.severity = "moderate"
        victimInjury.severityScore = 3
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "FMP"
        victimInjury.injuryLocation = "chest"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.6
    end
    if injuryNumber == 5
        victimInjury.name = "shrapnel blast leg"
        victimInjury.severity = "moderate"
        victimInjury.severityScore = 3
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "ED"
        victimInjury.injuryLocation = "lower extremity"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.2
        
    end
    if injuryNumber == 6
        victimInjury.name = "3rd degree burns > 50% body"
        victimInjury.severity = "moderate"
        victimInjury.severityScore = 3
        victimInjury.role =["role3"]
        victimInjury.treatmentLocation = "ED"
        victimInjury.injuryLocation = "external"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.1
    end
    if injuryNumber == 7
        victimInjury.name = "leg crush"
        victimInjury.severity = "moderate"
        victimInjury.severityScore = 3
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "ANY"
        victimInjury.injuryLocation = "lower extremity"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.1
    end
    if injuryNumber == 8
        victimInjury.name = "pneumothorax"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 4
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "ANY"
        victimInjury.injuryLocation = "chest"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 1.0
    end
    if injuryNumber == 9
        victimInjury.name = "penetrating shrapnel debris to head"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 4
        victimInjury.role =[""]
        victimInjury.treatmentLocation = "ED"
        victimInjury.injuryLocation = "head/neck"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.7
    end
    if injuryNumber == 10
        victimInjury.name = "bilateral leg amputation"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 5
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "FMP"
        victimInjury.injuryLocation = "lower extremities"
        victimInjury.hospitalNeed = "S"
        victimInjury.equipmentMitigation = 0.0
    end
    if injuryNumber == 11
        victimInjury.name = "smoke inhalation burns"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 4
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "ANY"
        victimInjury.injuryLocation = "chest"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.0 # 1.0 if CBRN protection
    end
    if injuryNumber == 12
        victimInjury.name = "Penetrating blast debris"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 4
        victimInjury.role =[""]
        victimInjury.treatmentLocation = "ED"
        victimInjury.injuryLocation = "chest"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.3
    end
    if injuryNumber == 13
        victimInjury.name = "arm amputation"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 4
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "FMP"
        victimInjury.injuryLocation = "upper extremity"
        victimInjury.hospitalNeed = "S"
        victimInjury.equipmentMitigation = 0.0
    end
    if injuryNumber == 14
        victimInjury.name = "spinal cord crush"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 4
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "ANY"
        victimInjury.injuryLocation = "spine"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.0
    end
    if injuryNumber == 15
        victimInjury.name = "Rupture of hollow viscera"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 4
        victimInjury.role =[""]
        victimInjury.treatmentLocation = "ED"
        victimInjury.injuryLocation = "abdomen"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.4
    end
    if injuryNumber == 16
        victimInjury.name = "unilateral leg amputation"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 4
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "FMP"
        victimInjury.injuryLocation = "upper extremity"
        victimInjury.hospitalNeed = "S"
        victimInjury.equipmentMitigation = 0.0
    end
    if injuryNumber == 17
        victimInjury.name = "asphyxia from smoke"
        victimInjury.severity = "moderate"
        victimInjury.severityScore = 3
        victimInjury.role = [""]
        victimInjury.treatmentLocation = "ANY"
        victimInjury.injuryLocation = "head/neck"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.0
    end
    if injuryNumber == 18
        victimInjury.name = "cervical fracture with transection"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 4
        victimInjury.role = ["role3"]
        victimInjury.treatmentLocation = "ANY"
        victimInjury.injuryLocation = "head/neck"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.0
    end
    if injuryNumber == 19
        victimInjury.name = "eye globe rupture"
        victimInjury.severity = "moderate"
        victimInjury.severityScore = 3
        victimInjury.role = ["role3"]
        victimInjury.treatmentLocation = "ANY"
        victimInjury.injuryLocation = "head/neck"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.0
    end
    if injuryNumber == 20
        victimInjury.name = "abdominal perforation"
        victimInjury.severity = "moderate"
        victimInjury.severityScore = 3
        victimInjury.role = ["role3"]
        victimInjury.treatmentLocation = "ANY"
        victimInjury.injuryLocation = "abdomen"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.3
    end
    if injuryNumber == 21
        victimInjury.name = "3rd degree burns > 75% body"
        victimInjury.severity = "severe"
        victimInjury.severityScore = 4
        victimInjury.role =["role3"]
        victimInjury.treatmentLocation = "ED"
        victimInjury.injuryLocation = "external"
        victimInjury.hospitalNeed = "ICU"
        victimInjury.equipmentMitigation = 0.0
    end


    return victimInjury
end

function clinconds(env :: Environment,patient :: Patient,SS_data :: Array,consts::Constants)
        
    if (consts.creator === true)
        d = patient.timeOfDeath
        ip = 7
        issLevel = patient.iss
        vicid = patient.vicid
        g=0.2
        if patient.isBleeding ==1
            g=1.0
        end
        
        #trauma part
       
        if issLevel >= 5 && issLevel <8
            g=0.5
        end
        if issLevel >= 8 && issLevel < 10
            g=0.6
        end
        if issLevel >= 10 && issLevel < 20
            g=0.8
        end        
        if issLevel >= 20  && issLevel < 26
            g=0.9
        end
        if issLevel >= 26
            g=1.0
        end
        
        
        
        
        if issLevel < 26

                if patient.age <= 12
                    b = 0.0055*patient.age -0.09
                end
                if patient.age > 70
                    b = -0.0067*patient.age +0.4433
                end
                if patient.age > 12 && patient.age <= 70
                    b = -0.03
                end
        else
            if patient.age <= 12
                b = 0.0055*patient.age -0.09 - 0.04
            end
            if patient.age > 70
                b = -0.0067*patient.age +0.4433 - 0.04
            end
            if patient.age > 12 && patient.age <= 70
                b = -0.03 
            end

        end
        if issLevel >= 8 
            a = b*d + 2.71828 #euler's number
            if a > 0
                a=a*(-1)
            end
        else 
            a = -8 
        end
        
       
        
        
        #a_rad,g_rad,radDeath,dose,dum2 = deltaR(env,0.0,patient,patient.distFromBlast,5000.0)
        patient.dose = 0.0
        g_rad = 1.0
        a_rad = 0.0
        radDeath = 1e6
        b_rad = -0.03
        SimParams = [a,b,g,d,ip,a_rad,b_rad,g_rad,radDeath]
        
        if consts.writeSimScore
            for t in 1:1000
           
                time_from_onset = convert(Float64,t)
                SimScore_all = generateSimedisScoreGompertz(env,time_from_onset*60,SimParams)
                # SimScore_phys = generateSimedisScoreGompertz(env,time_from_onset*60,SimParams,"Physical")
                # SimScore_chem = generateSimedisScoreGompertz(env,time_from_onset*60,SimParams,"Chemical")

                displayTime = time_from_onset
                SS_data = push!(SS_data,[vicid  displayTime  SimScore_all])
                t += 60
            end
        end
    else
        SimParams = patient.SimParams
       
        ip = patient.SimParams[5]
        vicid =patient.vicid
        if consts.writeSimScore
            for t in 1:1000
         
                time_from_onset = convert(Float64,t)
                SimScore = generateSimedisScoreGompertz(env,time_from_onset*60,SimParams)
                displayTime = time_from_onset
                SS_data = push!(SS_data,[vicid  displayTime  SimScore])
               
            end
        end
    # params = push!(params,[vicid  round(SimParams[1],digits =2) round(SimParams[2],digits=2) round(SimParams[3],digits=2)   round(SimParams[4],digits=2)  SimParams[5]])
        
    end
    return SimParams
end

function treatmentFunction(patient::Patient,t::Float64,tstart::Float64,tstop::Float64)
    #compute the deltaTreatment from the treatmentParams
   
    if t >= tstart
        
        deltaT = sqrt(patient.treatmentParams[1]*(t-tstart)/60.0) + exp(-patient.treatmentParams[2]*(t-tstart)/60.0)*sqrt(patient.treatmentParams[3]*(t-tstart)/60.0)
       

        return deltaT
    else 
        return 0.0
    end
        return 0.0
end

function distance(lat1, lon1, lat2, lon2)
    return sqrt((lat1 - lat2)^2 + (lon1 - lon2)^2)
end
# Function to find the closest lat, lon in the filtered dictionary
function find_closest_concentration(lat, lon, filtered_dict)
    closest_key = nothing
    min_dist = Inf
    closest_concentration = missing

    for (key, value) in filtered_dict
        key_lat, key_lon = key
        dist = distance(lat, lon, key_lat, key_lon)
        if dist < min_dist
            min_dist = dist
            closest_key = key
            closest_concentration = value
        end
    end

    return closest_concentration
end
function calculateDose(data_array::Array,filtered_data::Dict)
    # Step 8: Calculate cumulative dose in Gray
    global total_dose = 0.0
    dose = 0
    for d in data_array
        concentration = find_closest_concentration(d[1], d[2], filtered_data)
        if concentration !== missing
            # Convert concentration to dose using the conversion factor
            dose = 1e-6*(-3.02e-7*(concentration*concentration) +7.7e-3*concentration + 5e-2)/60.0 #data from Scherb and Hayashi for Fukushima µSv/h to mins then µSv to Gy = 1e-6
            
            global total_dose += dose
           
        end
    end

    return total_dose
end
function haversine_distance(start_coord::Tuple{Float64, Float64}, end_coord::Tuple{Float64, Float64})
    radius = 6371.0  # Earth's radius in kilometers
    
    lat1, lon1 = start_coord
    lat2, lon2 = end_coord
    
    # Convert degrees to radians
    φ1 = lat1 * π / 180.0
    φ2 = lat2 * π / 180.0
    Δφ = (lat2 - lat1) * π / 180.0
    Δλ = (lon2 - lon1) * π / 180.0
    
    # Haversine formula
    a = sin(Δφ/2)^2 + cos(φ1) * cos(φ2) * sin(Δλ/2)^2
    c = 2 * atan(sqrt(a) / sqrt(1 - a))
    
    return radius * c  # Distance in kilometers
end
function interpolate_path(start_coord::Tuple{Float64, Float64}, 
        end_coord::Tuple{Float64, Float64}, 
        speed_kmph::Float64, 
        start_time::Int64)

    # Calculate the distance between the start and end points
    total_distance_km = haversine_distance(start_coord, end_coord)

    # Calculate the total time required in minutes
    total_time_min = (total_distance_km / speed_kmph) * 60.0

    # Determine the number of points (every minute)
    num_points = Int(round(total_time_min))

    # Create latitude and longitude ranges
    latitudes = range(start_coord[1], stop=end_coord[1], length=num_points)
    longitudes = range(start_coord[2], stop=end_coord[2], length=num_points)

    # Generate timestamps in seconds
    timestamps = [start_time + (i-1) * 60 for i in 1:num_points]

    # Generate the path with latitude, longitude, and timestamp
    path = [(latitudes[i], longitudes[i], timestamps[i]) for i in 1:num_points]

    return path
end
function ParseConc(json_file::String)
  
    json_data = JSON.parsefile(json_file)
    # Step 2: Convert string keys to tuples and build a dictionary

    filtered_data = Dict{Tuple{Float64, Float64}, Float64}()

    # Step 3: Parse and filter the data for time = 259200
    for entry in json_data
        for (key_str, value) in entry
            # Parse the string key into a tuple
            tuple_key = eval(Meta.parse(key_str))  # Converts the string key to a tuple
            lat, lon = tuple_key  # Unpack the tuple into lat, lon
            
            # Debugging: Print the parsed values
            filtered_data[(lat, lon)] = value
            
       
        end
    end
    return filtered_data
end

function getLD50(condition::Int)
    if condition == 0
        return 4.5#in Gy no med treatment
    end
    if condition ==1 #medical treatment ok except G-CSF
        return 6.8
    end
    if condition ==2 #granulocyte colony-stimulating factor treatment
        return 8.5
    end
end

function generateTraj(origin:: OpenStreetMapX.LLA,destination::OpenStreetMapX.LLA,transport::String,m) #genere aleatoirement traj pour patient partant de origin
    
    traj=[]
    route=[]
    originNode = point_to_nodes((origin.lat,origin.lon),m)
    destinationNode = point_to_nodes((destination.lat,destination.lon),m)
    route,_,_= OpenStreetMapX.shortest_route(m,originNode,destinationNode)
    for k in 1:length(route)
        locs = [OpenStreetMapX.LLA(m.nodes[n],m.bounds) for n in route[k]]
        push!(traj, [(loc.lat, loc.lon, k) for loc in locs ])
    end

    return traj
      

end


function whole_body_radiation(distance::Float64, yield::Float64) #this is for a tac nuke explosion lookup flexpart instead
    tnt_energy = 4.184e9 # Energy released per ton of TNT in joules
    joules_per_gray = 1e-2 # Joules of energy required to produce 1 Gray of whole-body radiation
    energy = yield * tnt_energy
    radius = 10 * (energy / joules_per_gray)^(1/3) # Blast radius in meters
    dose_rate = energy / (4 * π * distance^2) # Dose rate in joules per second per square meter
    dose=dose_rate * (1 - exp(-7.8e-3 * distance^(1.08 - 0.075 * log10(distance)) * radius^(-0.82))) # Total whole-body radiation in Grays
    return dose
end
function plotrad(yield::Float64)
   #plots the release vs distance for a yield in Tonnes of TNT
    # Define the range and yield
    distance_range = 0:25.0:100000.0 #up to 100km


    # Calculate the radiation dose at each distance
    dose = [whole_body_radiation(d, yield) for d in distance_range]

    # Filter the distances with doses below 8.5 Gy
    #dose_range = collect(filter(d -> Simedis.whole_body_radiation(d, yield) < 8.5, distance_range))

    # Create a plot of the dose versus distance
    # plot(distance_range, dose, xlabel = "Distance (m)", ylabel = "Dose (Gy)", label = "Radiation Dose")
    # scatter!(dose_range, [Simedis.whole_body_radiation(d, yield) for d in dose_range], label = "Dose below 8.5 Gy")
    # savefig("rad_plot.png")
    
    output_file = "dose_vs_distance.txt"
    data = [distance_range dose]
    writedlm(output_file, data)
end
function deltaR(env:: Environment, t::Float64, patient :: Patient,distance::Float64,yield::Float64)
    #duration in hours, dose in Gy, LD50 in Gy function of treatment
    #yield in Joules
    dose = whole_body_radiation(distance::Float64, yield::Float64) 
    duration=3.0
    LD50 = getLD50(0)
    D_death = LD50 /(-0.2351*0.8946^(dose/duration)*(dose/duration)^-0.2876 + 0.9947)
    t_death = 429*dose^(-1.3) #in hours
    if dose >= D_death && t >= t_death*3600
       
        gamma=1.0
        b=-0.03
        a_new = b*(t_death*60.0) + 2.7183  #t in mins
        dR = 20.0 #equivalent to 75% of max SIMSCORE and incapacitated following AMedP-8
        
        return a_new,gamma,t_death*60,D_death,dR
    end
    if dose >= D_death && t < t_death*3600
        gamma=0.9  
        
        b=-0.03     
        a_new = b*(t_death*60.0) + 2.7183  #t in mins
        dR = 20.0 - (20.0 - (20.0-20.0*(exp(-exp(a_new-0.03*t/60.0))))^gamma) #values 20 if gompertzRad = 0, 0 if gompertz rad = 20, t in seconds conv to mins
        
        return a_new,gamma,t_death*60,D_death,dR
    #takes the dose in Gy and converts it into a delta(SS)
    
    end
    if dose < D_death && t < t_death*3600
        gamma=0.8
        b=-0.03
        a_new = b*(t_death*60.0) + 2.7183  #t in mins
        dR = 20.0 - (20.0 - (20.0-20.0*(exp(-exp(a_new-0.03*t/60.0))))^gamma) #values 20 if gompertzRad = 0, 0 if gompertz rad = 20, t in seconds conv to mins
        
        return a_new,gamma,t_death*60,D_death,dR
    end
    if dose < D_death && t >= t_death*3600
        gamma=1.0
        patient.SimParams[8] = 1.0
        b=-0.03
        a_new = b*(t_death*60.0) + 2.7183  #t in mins
        dR = 20.0 - (20.0 - (20.0-20.0*(exp(-exp(a_new-0.03*t/60.0))))^gamma) #values 20 if gompertzRad = 0, 0 if gompertz rad = 20, t in seconds conv to mins
        
        return a_new,gamma,t_death*60,D_death,dR
    end

end
function deltaIP(ipLevel::Int64,t::Float64,consts::Constants)

    if consts.tabularIP  
        if ipLevel == 1 || ipLevel == 2 || ipLevel == 7
                
                return 0.0 #IP1&2 have no noticeable symptoms
                
        end
        if ipLevel == 3
            if t<=3
                return 0.0 
            else
                return 1.0
            end
        end
        ip4=[0.0,1.0,6.0,3.0,2.0,1.0,0.0]
        tx=[1.0,5.0,7.0,21.0,106.0,1006.0,21000.0]
        ip5=[0.0,1.0,6.0,8.0,6.0,3.0,2.0,1.0,0.0]
        tx1=[1.0,5.0,7.0,13.0,23.0,68.0,248.0,1008.0,20168.0]
        ip6=[0.0,1.0,6.0,15.0,20.0,20.0,20.0,20.0]
        tx2=[1.0,5.0,7.0,11.0,15.0,25.0,70.0,20000.0]

        if ipLevel == 4
            for i in 1:length(tx)-1
                if t >= tx[i] && t < tx[i+1]
                   
                    return ip4[i] 
                
                end
            end
            return 0.0
        end
        if ipLevel == 5
            for i in 1:length(tx1)-1
                if t >= tx1[i] && t < tx1[i+1]
                    return ip5[i] 
                
                end
            end
            return 0.0
        end
        if ipLevel == 6
            for i in 1:length(tx1)-1
                if t >= tx2[2] && t < tx2[i+1]
                    return ip6[i] 
            
                end
            end
            return 0.0
        end
    else

        if ipLevel == 1 || ipLevel == 2 || ipLevel == 7
                
            return 0.0 #IP1&2 have no noticeable symptoms
            
        end
        if ipLevel == 3
            if t<5
                return 0.0 
            else
                return 1.0
            end
        end
        if ipLevel == 4
            if t<5
                return 0.0 #asymptotic behavior and no symptoms before 5 mins anyway discard chem effects
            else
                if t==5
                    return 1.0
                else
                    if t > 1000.0
                        return 0.0

                    else
                        return (3.35 +(10.2/(2^(2.17/2.0)*gamma(2.17/2.0))*t^((2.17/2.0)-1)*exp(-t/180.0))^1.0)
                    end
                end
            end
        end
        if ipLevel == 5
            if t<5
                return 0.0 #asymptotic behavior and no symptoms before 5 mins anyway discard chem effects
            else
                return 1.0 + (21/(2^(1.45/2.0)*gamma(1.45/2.0))*t^((1.45/2.0)-1)*exp(-t/80.0))^0.8
            end
        end
        if ipLevel == 6
            if t<5
                return 0.0 #asymptotic behavior and no symptoms before 5 mins anyway discard chem effects
            else
             return 20*exp(-exp(4-0.48*t))
            end
        end
       
    end

end
function checkHealth(env :: Environment, patient :: Patient,t:: Float64,thresholdHealth::Float64)

    
    SimParams = patient.SimParams
    y=generateSimedisScoreGompertz(env,t*60,SimParams)

    if y <= thresholdHealth
       
        return 1
    else
        return 0

    end
end
function generateSimedisScoreGompertz(env:: Environment,t:: Float64, SimParams :: Array, subType::String="Both" )
    #input t is in seconds!
        b = SimParams[1]
        c= SimParams[2]
        gamma = SimParams[3]
        ip = SimParams[5]
        f = 0.0
   
        return y = 20.0 - (20.0-20.0*(exp(-exp(b-c*t/60.0))))^gamma+f*sin(((rand(Float64)))*pi/2-pi/2)
       
end
function generateSimedisScoreGompertzPatient(env:: Environment,t:: Float64, patient::Patient, consts::Constants,subType::String="Both")
    #input t is in seconds!
        b = patient.SimParams[1]
        c = patient.SimParams[2]
        gamma = patient.SimParams[3]
        ip = patient.SimParams[5]
        maxSimScore = patient.maxSimScore 
           
        ipI = convert(Int64,ip) 
        deltaC =deltaIP(ipI,t/60.0+now(env),consts)
        dR1=0.0
        #a,a1,a2,a3,dR1=deltaR(env,t/60.0,patient,patient.distFromBlast,5000.0)
        dT = treatmentFunction(patient,t/60.0,patient.treatmentParams[4]/60.0,patient.treatmentParams[5]/60.0)
      
        
        
        f=0.0
        
        if subType == "Physical"
            deltaC =0
            dR1=0
       
        end
        if subType == "Nuclear"
            deltaC =0
            y = maxSimScore - dR1
        end
        if subType == "Chemical"
            dR1=0
            y = maxSimScore - deltaC
        end
        if patient.isBleeding == 1
            bv=computeBloodVolume(patient,t/3600.0,0.0,0.5,0.08,0.01,false)
            patient.bloodvolume=bv
            if bv<0.2
                return 0.0 #patient died by hemorrhage
            end
        end
        if patient.isBleeding == 0
            bv=computeBloodVolume(patient,t/3600.0,0.0,2.0,0.03,0.01,true)
            patient.bloodvolume=bv
            if bv<2.0
                gamma=1.0
            end
        end
                
        y = max(0.0,min(20.0,maxSimScore - (maxSimScore-maxSimScore*(exp(-exp(b-c*((t/60.0))))))^gamma - deltaC -dR1+dT))
       
           
       
        if y < 17
            patient.incapacitated =0 #remove references to incapacitated for now
        else
            patient.incapacitated =0
        end
        
        
        return y
end

function applyTreatment(env :: Environment, patient :: Patient,t:: Float64,treatmentTime :: Float64,treatmentType::Int64,treatmentStartTime::Float64,treatment_data::Array,paramset::Parameterset,consts)
    #t is seconds!but generatesscore takes seconds and returns a t in minutes!
   
    patient.treatmenthasstarted = true
    simT = generateSimedisScoreGompertzPatient(env,t,patient,consts)
    #write a function for the treatment params (function of the ISS)
    if patient.iss > 0  
        deltaSS = 20.0
    else
        deltaSS = 20.0 #look up/discuss
    end 
    

    # deltaSS1=0.3*deltaSS #fmp
    deltaSS3 = 1.0*deltaSS #hospital

   
    # if treatmentType == 1 
            
    #         if t/60.0 == treatmentTime/60.0 + treatmentStartTime#treatmentFinished adapting the parameters
    #             SimParams = patient.SimParams
    #             b = SimParams[1]
    #             c= SimParams[2] 
    #         end
    #         if t/60.0 == treatmentTime/60.0 + treatmentStartTime#treatmentFinished adapting the parameters
    #             SimParams = patient.SimParams
    #             b = SimParams[1]
    #             c= SimParams[2] 
                
    #             c= c + 0.02 #fmp C increase
    #             SimParamsPostTreatment = [b,c,SimParams[3],SimParams[4],SimParams[5],SimParams[6],SimParams[7],SimParams[8],SimParams[9]]
    #             patient.SimParams = SimParamsPostTreatment
    #             patient.treatmentEndTime += t/60.0
    #             if patient.isBleeding == 1
    #                 patient.isBleeding = 0
    #                 patient.SimParams[3]=0.9
    #             end
                
    #             return 1.0
            
    #         end
   
    # end
    
  
    # if treatmentType == 2
      
        
    #     if t/60.0 >= (treatmentTime/60.0 + treatmentStartTime)
      
            
    #         if t/60.0 == (treatmentTime/60.0 + treatmentStartTime)
           
                
    #             SimParams = patient.SimParams
    #             b = SimParams[1]
    #             c= SimParams[2] 
    #             #bPost = b - 2.0
    #             if paramset.TranspsupervisionLevel == "Normal"
    #             SimParams = patient.SimParams
    #             b = SimParams[1]
    #             c= SimParams[2] 
    #             #bPost = b - 2.0
    #             if paramset.TranspsupervisionLevel == "Normal"
                    
    #                 c= c + 0.02
    #                 c= c + 0.02
                
    #             end
    #             if paramset.TranspsupervisionLevel == "Low"
    #                 c= c + 0.01
               
    #             end
    #             if paramset.TranspsupervisionLevel == "Low"
    #                 c= c + 0.01
    #             end

    #             SimParamsPostTreatment = [b,c,SimParams[3],SimParams[4],SimParams[5],SimParams[6],SimParams[7],SimParams[8],SimParams[9]]
    #             patient.SimParams = SimParamsPostTreatment
    #             patient.treatmentEndTime += t/60.0
             
    #             return 1.0
              
    #         end
    #         if t/60.0 > (treatmentTime/60.0 + treatmentStartTime)
    #             simT = generateSimedisScoreGompertzPatient(env,t-treatmentStartTime*60,patient,consts)
    #             return simT
    #         end
           
           

    #     end
    
    # end

    if treatmentType == 3 #lifeSavingTreatment
          
          
            if t/60.0 >= (treatmentTime/60.0 + treatmentStartTime)
                
                if t/60.0 == (treatmentTime/60.0 + treatmentStartTime)
                    patient.SimParams[3]=0.2
                    y = generateSimedisScoreGompertzPatient(env,t,patient,consts)
                 
                 
                    simT =  min(20.0,y+deltaSS3)
                    patient.maxSimScore = simT
                    return simT

                end
                if t/60.0 > (treatmentTime/60.0 + treatmentStartTime)
                    
                    
                    patient.SimParams[3]=0.2
                    y = generateSimedisScoreGompertzPatient(env,t,patient,consts)
                    simT =  min(20.0,y+deltaSS3)
               
                    patient.isLethallyInjured = 0  #saved
                    return simT
                end

            end
    end

end
function setTraumaTriage(env::Environment,patient::Patient,consts)
    
   
    TraumaScore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
   
    
    #iss = patient.iss
   # canwalk = patient.mobile
    if patient.triage == 9
        return 9
    end
  
    if  checkMobility(env,patient) && (TraumaScore >= 17.0) && patient.isBleeding == 0
        
        patient.triage = 2
        return 2
    elseif (TraumaScore < 17.0)

        patient.triage = 1
        return 1

    else
        patient.triage = 2
        return 2
            
    end
    
end
function setRadioTriage(env::Environment,patient::Patient,consts::Constants)
    RadioScore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts,"Nuclear")
    #iss = patient.iss
   # canwalk = patient.mobile
    if patient.triage == 9
        return 9
    end
  
    if  checkMobility(env,patient) && (RadioScore >= 17.0)
        
        patient.triage = 3
        return 3
    elseif (RadioScore < 17.0)

        patient.triage = 1
        return 1

    else
        patient.triage = 2
        return 2
            
    end


end
function setTriage(env::Environment,patient::Patient,consts::Constants)
  
    
    triage = min(setTraumaTriage(env,patient,consts),5,setRadioTriage(env,patient,consts))
    if patient.triage !=9
        patient.triage = triage
    else #if dead stay dead
        patient.triage = 9
        patient.isdead = 1
    
    end

    return triage

end
function triageMistakes(triage :: Int64, mistake :: String)
    # returns triage category according to preset error: over (higher priority), under (lower priority), or correct
    if mistake == "over"
        if triage == 2 || triage == 3
            return (triage - 1)
        elseif triage == 4
            return 1
        else
            return triage
        end
    elseif mistake == "under"
        if triage == 1 || triage == 2
            return (triage + 1)
        else
            return triage
        end
    elseif mistake == "correct"
        return triage
    end
end
function DetermineTriageCorrectness(TriageIncorrectProb :: Tuple{Int64,Int64})

    probUnder = TriageIncorrectProb[1]
    probOver = TriageIncorrectProb[2]

    randomnr = rand()
    if randomnr < (probUnder/100)
        return "under"
    elseif randomnr > (1 - probOver/100)
        return "over"
    else
        return "correct"
    end
end
function isdead(patient ::Patient)
   if patient.triage == 9
    return true
   else
    return false

   end

end
function vicPriorityScore(env :: Environment,patient ::Patient,consts)
    triagecat = patient.triage

 
    simedisscore = trunc(Int64,generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts))
    
    return triagecat*100+simedisscore
end
function getTriage(issvalue::Int64,iplevel::Int64)
    triageISS = 5
    triageIP = 5
    if issvalue >= 25
        triageISS = 1
    end
    if issvalue < 5 && issvalue >= 3
        triageISS = 3
    end
    if issvalue >= 5 && issvalue < 25

        triageISS = 2
    end
    if issvalue < 1
        triageISS = 5
    end

    if iplevel >= 5
        triageIP = 1
    end
    if iplevel <= 2
        triageIP = 3
    end
    if iplevel > 6 #ip7 is 0.0 dose
        triageIP = 5
    end
    if iplevel < 5 && iplevel > 2
        triageIP = 2
    end

    triageVal = min(triageIP,triageISS)

    return triageVal
end
@resumable function releaseResourceDirect(env :: Environment, medic :: Resource, request :: SimJulia.Put)
    @yield request
    @yield release(medic)
end
"""
    scheduleMTFPersonnel(env:: Environment, mediclist::MedicList,paramset::Parameterset, inputdb::SQLite.DB,consts::Constants)

Mans each MTF with doctors and nurses to be used as resources for treatments
"""

@resumable function scheduleMTFPersonnel(env:: Environment, mediclist::MedicList,paramset::Parameterset, inputDB::SQLite.DB,consts::Constants)

    hospList = DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Hospital" ))
    for i in 1:nrow(hospList)
        if hospList[i,:].type == "R1"
            @process addMedic(env,mediclist.fmpDoctor,1)
            @process addMedic(env,mediclist.fmpNurse,4)
        end
        if hospList[i,:].type == "R2" 
            @process addMedic(env,mediclist.fmpDoctor,2)
            @process addMedic(env,mediclist.fmpNurse,6)
        end
        if hospList[i,:].type == "R2E"
            @process addMedic(env,mediclist.fmpDoctor,4)
            @process addMedic(env,mediclist.fmpNurse,12)
        end
        if hospList[i,:].type == "R2F"
            @process addMedic(env,mediclist.fmpDoctor,2)
            @process addMedic(env,mediclist.fmpNurse,6)
        end
        if hospList[i,:].type == "R2B"
            @process addMedic(env,mediclist.fmpDoctor,3)
            @process addMedic(env,mediclist.fmpNurse,8)

        end
        if hospList[i,:].type == "DECON"
            @process addMedic(env,mediclist.fmpDoctor,1)
            @process addMedic(env,mediclist.fmpNurse,6)
        end

    end

end
function timerandomised(time :: Real, label :: String,consts::Constants)
    if !consts.stochasTimes || time == 0.0
        return time
    elseif label == "Treatment"
        distributionName = consts.TreatTimeDist
        boundary = consts.TreatTimeBound
    elseif label == "Travel"
        distributionName = consts.TravTimeDist
        boundary = consts.TravTimeBound
    elseif label == "Rescue"
        distributionName = consts.RescueTimeDist
        boundary = consts.RescueTimeBound
    else
        distributionName = consts.OtherTimeDist
        boundary = consts.OtherTimeBound
    end

    if distributionName == "Normal"
        sig = boundary*time/consts.sigmaN
        d = Truncated(Normal(time,sig), time*(1-boundary),time*(1+boundary))
        return rand(d)
    elseif distributionName == "LogNormal"
        mu = log(boundary*time)
        d = Truncated(LogNormal(mu,sqrt(consts.sigmaLN)),0,time*(1+boundary))
        return (rand(d) + time * (1 - boundary))
    elseif distributionName == "Uniform"
        d = Uniform(time*(1-boundary),time*(1+boundary))
        return rand(d)
    elseif distributionName == "Triangular"
        d = TriangularDist(time*(1-boundary),time*(1+boundary),time)
        return rand(d)
    else
        return time # to catch not defined distributions
        println("ERROR IN timerandomised: DISTRIBUTION NOT DEFINED, label = $label")
    end
end
function createHospital(hospitalProps::Array,locationx::Float64,locationy::Float64,m)
    hospNode = point_to_nodes((hospitalProps[16],hospitalProps[17]),m)
    originNode = point_to_nodes((locationx,locationy),m)
    route, distance,Time = OpenStreetMapX.fastest_route(m,hospNode,originNode)            
    driveTime = (Time/60.0)*1.7
    
    dfHospital = DataFrame(ID = hospitalProps[1],Name= hospitalProps[2], CapT1_A_Medium=hospitalProps[3],CapT2_A_Medium=hospitalProps[4],CapT3_A_Medium=hospitalProps[5],CapT1_P_Medium=hospitalProps[6],CapT2_P_Medium=hospitalProps[7],CapT3_P_Medium=hospitalProps[8],DriveTime=driveTime,SurgeCapT1_A=hospitalProps[10],SurgeCapT2_A=hospitalProps[11],SurgeCapT3_A=hospitalProps[12],SurgeCapT1_P=hospitalProps[13],SurgeCapT2_P=hospitalProps[14],SurgeCapT3_P=hospitalProps[15])
    newHospital = Hospital(dfHospital[1,:],"Medium")
   
    return newHospital
end
function defineHospitals(hospitalProps::Array,locationx::Float64,locationy::Float64,m,hospitals::DataFrame)
    hospitalQueue = Array{Hospital}(undef,size(hospitals,1))
    for i = 1:size(hospitals,1)
        hospNode = point_to_nodes((hospitals[i,:].LAT,hospitals[i,:].LON),m)
        originNode = point_to_nodes((locationx,locationy),m)
        route, distance,Time = OpenStreetMapX.fastest_route(m,hospNode,originNode)  
        driveTime = (Time/60.0)*1.7
        dfHospital = hospitals[i,:]
        hospital = Hospital(dfHospital,"Medium")
        hospitalQueue[i] = hospital
        
    end
    return hospitalQueue
end

function CreateHospitalQueue(inputDB :: SQLite.DB, locations,m,routes,paramset,offRoad::Bool) #hospitalQueue per threat

    hospitals = DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Hospital" ))
    
    
    hospitalQueue = Array{Hospital}(undef,size(DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Hospital" )),1))
    
    driveTime=[0.0, 0.0 ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    for i = 1:size(hospitals,1)
        dfHospital = hospitals[i,:]
        hospNode = point_to_nodes((hospitals[i,:].LAT,hospitals[i,:].LON),m)
        for k in 1:length(locations)
            originNode = point_to_nodes((locations[k].lat,locations[k].lon),m)
            route, distance,Time = OpenStreetMapX.shortest_route(m,hospNode,originNode) 
            if offRoad
                route_, distance_,Time_ = OpenStreetMapX.fastest_route(m,hospNode,originNode) 
              
               # push!(routes,route_)
                originECF = ECEFfromLLA(wgs84)(Geodesy.LLA(locations[k].lat,locations[k].lon))
                hospECF = ECEFfromLLA(wgs84)(Geodesy.LLA(hospitals[i,:].LAT,hospitals[i,:].LON))
            
                distToHosp=OpenStreetMapX.distance(originECF.x,originECF.y,originECF.z,hospECF.x,hospECF.y,hospECF.z)
               
                # if paramset.simurun ==1
                #     push!(routes,route)
                # end OFFROAD -> no routes
                driveTime[k] = (distToHosp*2.0/(30*(16.6666))) #time in mins is distance in m divided by speed in m/min
            
                dfHospital.DriveTime = driveTime[k]
            else
                route_, distance_,Time_ = OpenStreetMapX.fastest_route(m,hospNode,originNode) 
                if paramset.simurun ==1
                    push!(routes,route)
                end
                driveTime[k] = (Time/60.0)*1.7
                dfHospital.DriveTime = driveTime[k]

            end
        end
       
        for k in 1:length(locations)
            originNode = point_to_nodes((locations[k].lat,locations[k].lon),m)
            route, distance,Time = OpenStreetMapX.shortest_route(m,hospNode,originNode) 
            if offRoad
                route_, distance_,Time_ = OpenStreetMapX.fastest_route(m,hospNode,originNode) 
              
               # push!(routes,route_)
                originECF = ECEFfromLLA(wgs84)(Geodesy.LLA(locations[k].lat,locations[k].lon))
                hospECF = ECEFfromLLA(wgs84)(Geodesy.LLA(hospitals[i,:].LAT,hospitals[i,:].LON))
            
                distToHosp=OpenStreetMapX.distance(originECF.x,originECF.y,originECF.z,hospECF.x,hospECF.y,hospECF.z)
               
                # if paramset.simurun ==1
                #     push!(routes,route)
                # end OFFROAD -> no routes
                driveTime[k] = (distToHosp*2.0/(30*(16.6666))) #time in mins is distance in m divided by speed in m/min
            
                dfHospital.DriveTime = driveTime[k]
            else
                route_, distance_,Time_ = OpenStreetMapX.fastest_route(m,hospNode,originNode) 
                if paramset.simurun ==1
                    push!(routes,route)
                end
                driveTime[k] = (Time/60.0)*1.7
                dfHospital.DriveTime = driveTime[k]

            end
        end
       
        hospital = Hospital(dfHospital,"Medium")
        hospitalQueue[i] = hospital
    end
 
    return routes,hospitalQueue
 
    return routes,hospitalQueue
end
function findHospitalswithFacility(inputDB ::SQLite.DB, facility ::String, adult ::Bool)
    if facility == "H"
        facility = "OPC"
    end
   
    if adult 
    
        df=DataFrame(DBInterface.execute( inputDB, "SELECT ID FROM Hospital WHERE "* facility *"_A = 1"))[!,1]
        return df
    else
    
        df=DataFrame(DBInterface.execute( inputDB, "SELECT ID FROM Hospital WHERE "* facility *"_P = 1" ))[!,1]
        return df
        
    end
end
function findHospitalTravelTimes(hospitalQueue ::Array{Hospital,1}, inputDB ::SQLite.DB, hospitalneed ::String, adult ::Bool,locationx::Float64,locationy::Float64,m,quadrant,offroad::Bool)

   
    if hospitalneed == "H"
        hospitalneed = "OPC"
    end
    if adult
        facility = hospitalneed*"_A"
        
    else
        facility = hospitalneed*"_P"
    end
    travelTimes=[zeros(Float64,length(hospitalQueue)) zeros(length(hospitalQueue))]
    if length(hospitalQueue) > 0
    
            
           
        result = DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Hospital WHERE $facility = 1" )) #hospital who can admit victim
        lla_location = LLAfromUTM(quadrant,true,wgs84)
        location_lla = lla_location(Geodesy.UTM(locationx,locationy,0.0))

        if !offroad
               
                
            for i in 1:size(result,1)
                
                hospNode = point_to_nodes((result[i,:].LAT,result[i,:].LON),m)
                originNode = point_to_nodes((location_lla.lat,location_lla.lon),m)
                route, distance,Time = OpenStreetMapX.fastest_route(m,originNode,hospNode) 
                driveTime = (Time/60.0)*1 #to calibrate with traffic info
                #result[i,:].DriveTime = driveTime
                travelTimes[i,:] = [driveTime convert(Int64,result[i,:].ID)]
              
            
            end
        else
            for i in 1:size(result,1)
                originECF = ECEFfromLLA(wgs84)(Geodesy.LLA(location_lla.lat,location_lla.lon))
                hospECF = ECEFfromLLA(wgs84)(Geodesy.LLA(result[i,:].LAT,result[i,:].LON))
            
                distToHosp=OpenStreetMapX.distance(originECF.x,originECF.y,originECF.z,hospECF.x,hospECF.y,hospECF.z)
            #* println("MTF: $(hospitals[i,:].Name) dist: $distToHosp drivetime: $((distToHosp*1.7/(30*(16.6666)))) from location $k")
                # if paramset.simurun ==1
                #     push!(routes,route)
                # end OFFROAD -> no routes
                driveTime = (distToHosp*1.0/(50*(16.6666))) #time in mins is distance in m divided by speed in m/min
                
                travelTimes[i,:] = [driveTime convert(Int64,result[i,:].ID)]
            
            end

        end

            return travelTimes
    
    else 
        return missing
    end
    
    return 
end
@resumable function scheduleAll(env :: Environment, firefighters :: Resource, mediclist :: MedicList, boatlist:: BoatList,ambulist :: AmbuList, medevac::Medevac,transportqueue ::Array,transportqueueT3 ::Array, hospitalQueue ::Array{Hospital,1}, paramset :: Parameterset, inputDB :: SQLite.DB,consts::Constants,treatment_data,outputDB,disasters,m,evac::Bool,offroad::Bool,homeUTM)
    if evac
        @process scheduleResetSimulator(env, hospitalQueue, 9999)
        # @process scheduleMMT(env, mediclist, ambulist, paramset, inputDB,consts,outputDB)
        #@process scheduleMEDEVAC(env, mediclist, medevac, paramset, inputDB,consts,outputDB)
        @process scheduleBOAT(env, boatlist, paramset, inputDB,consts,outputDB,disasters,m,0)
   
       #Æ @process boatTransport(env ::Environment, transportqueue ::Array, boatlist ::BoatList, paramset ::Parameterset,consts::Constants,inputDB,outputDB,false,1)
    # @process checkTransport(env, transportqueue, transportqueueT3, mediclist, ambulist, hospitalQueue, paramset,treatment_data,consts,inputDB,outputDB,disasters,m)
    else
        @process scheduleResetSimulator(env, hospitalQueue, 9999)
        
        if length(disasters) >= 2
           # @process scheduleSR(env, firefighters, paramset,consts,outputDB,inputDB,disasters,m,0,offroad)
           # @process scheduleMMT(env, mediclist, ambulist, paramset, inputDB,consts,outputDB)
          #  @process scheduleMEDEVAC(env, mediclist, medevac, paramset, inputDB,consts,outputDB)
            @process scheduleAMB(env, ambulist,mediclist, paramset, inputDB,consts,outputDB,disasters,m,0,[60.0,500.0],offroad)
            @process scheduleMTFPersonnel(env, mediclist,paramset,inputDB,consts) #man mtfs with personnel
        
        end
        if length(disasters) == 1
            @process scheduleSR(env, firefighters, paramset,consts,outputDB,inputDB,disasters,m,0,0.0,offroad)
         #   @process scheduleMMT(env, mediclist, ambulist, paramset, inputDB,consts,outputDB)
           # @process scheduleMEDEVAC(env, mediclist, medevac, paramset, inputDB,consts,outputDB)
            @process scheduleAMB(env, ambulist, mediclist,paramset, inputDB,consts,outputDB,disasters,m,0,[0.0],offroad)

        end
        #@process boatTransport(env ::Environment, transportqueue ::Array, boatlist ::BoatList, paramset ::Parameterset,consts::Constants,inputDB,outputDB)
        @process checkTransport(env, transportqueue, transportqueueT3, mediclist, ambulist, hospitalQueue, paramset,treatment_data,consts,inputDB,outputDB,disasters,m,offroad,homeUTM)
    end
end
@resumable function scheduleResetSimulator(env ::Environment, hospitalQueue ::Array{Hospital,1}, timeouttime)
    @yield timeout(env, timeouttime)
    stopSurge(hospitalQueue)
end
@resumable function startBuildingFMP(env :: Environment, mediclist :: MedicList,consts)

    @yield timeout(env, max(timerandomised(consts.setupFMPtime, "Other",consts),0))
    # time for medics to setup FMP ! 0 if FMP is a seized Hosp
    mediclist.TreatmentAtFMP = true
end
@resumable function releaseMedicDirect(env :: Environment, medic :: Resource, request :: SimJulia.Put)
    @yield request
    @yield release(medic)
end
@resumable function addMedic(env :: Environment, medic :: Resource, number :: Int64)
    medic.capacity += number
    @yield request(medic, priority = 1)
    @yield release(medic)
end
@resumable function informresourceLastVictim(env ::Environment, resource ::Resource)
    @yield request(resource, priority = 9999)
    resource.capacity -= 1
    resource.level -= 1
end
@resumable function formmmt(env ::Environment, mmt ::Resource, doctor ::Resource, nurse ::Resource)
    if doctor.level < doctor.capacity && nurse.level < nurse.capacity
        # take one dr and one nu (no smaller level, dr and nurse are idle, capacity > 0, since min(level) = 0)
        doctor.capacity -= 1
        nurse.capacity -= 1
    end
    @process addMedic(env, mmt, 1)
end
@resumable function checkTransport(env ::Environment, transportqueue ::Array, transportqueueT3 ::Array, mediclist ::MedicList, ambulist ::AmbuList, hospitalQueue ::Array{Hospital,1}, paramset ::Parameterset,treatment_data::Array,consts::Constants,inputDB,outputDB,disasters,m,offroad,homeUTM)
    # triggered by release of ambulist.transportchecktrigger, then activated transport function to transport if possible, then calls itself again, waiting for next trigger
    @yield request(ambulist.transportchecktrigger)
    locationUTM =  Geodesy.UTM(disasters[1][1],disasters[1][2],0.0)
    locationx =locationUTM.x #start point for transport Lat
    locationy = locationUTM.y #start point for transport Lon
    
    @process transport(env, transportqueue, transportqueueT3, mediclist, ambulist, hospitalQueue, paramset,treatment_data,consts,inputDB,outputDB,locationx,locationy,m,offroad,homeUTM)
    @process checkTransport(env, transportqueue, transportqueueT3, mediclist, ambulist, hospitalQueue, paramset,treatment_data,consts,inputDB,outputDB,disasters,m,offroad,homeUTM)
end
function checkTotalCapacity(hospitalQueue ::Array{Hospital,1}, triagecat ::Int64, facility ::String, adult ::Bool,inputDB)
    # return added capacities of all hospitals for specified triagecat, facility, and agecategory (adult or peadiatric)
   
    if facility == "H"
        facility = "OPC"
    end
    totalcapacity = 0
    hospitalswithfacility = findHospitalswithFacility(inputDB, facility, adult)
    
    for hospital in hospitalQueue
        
        if adult && hospital.ID in hospitalswithfacility
            if triagecat == 1 || triagecat == 4
                totalcapacity += hospital.CapT1
            elseif triagecat == 2
                totalcapacity += hospital.CapT2
            elseif triagecat == 3
                totalcapacity += hospital.CapT3
            end
        elseif !adult && hospital.ID in hospitalswithfacility
            if triagecat == 1 || triagecat == 4
                totalcapacity += hospital.CapT1Ped
            elseif triagecat == 2
                totalcapacity += hospital.CapT2Ped
            elseif triagecat == 3
                totalcapacity += hospital.CapT3Ped
            end
        end

    end
  
    
    return totalcapacity
end
function checkT3capacity(hospitalQueue ::Array{Hospital,1}, hospitalID ::Int64, adult ::Bool)
    # checks T3 capacity
    for hospital in hospitalQueue
        if hospital.ID == hospitalID
            if adult && hospital.CapT3 > 0
                return true
            elseif !adult && hospital.CapT3Ped > 0
                return true
            else
                return false
            end
        end
    end
    println("no hospital with id $hospitalID")
end

function findClosestHospitalId(hospitalQueue ::Array{Hospital,1}, triagecat ::Int64, hospitalneed ::String, adult ::Bool,inputDB,locationx,locationy,m,quadrant,offroad)
    # returns ID and TravelTime of closest hospital for specified triagecat, facility, and agecategory (adult or peadiatric)
   
    traveltimes = findHospitalTravelTimes(hospitalQueue, inputDB, hospitalneed, adult,locationx,locationy,m,quadrant,offroad)
   
   
    sortedTT=sortslices(traveltimes,dims=1,by=x->x[1],rev=false)
  
      
    if length(traveltimes) > 0
        for i in 1:length(traveltimes)
            for hospital in hospitalQueue
                
                if hospital.ID == convert(Int32,(sortedTT[i,2]))
               
                  
                  
                    if adult
                        if (triagecat == 1 || triagecat == 4) && hospital.CapT1 > 0
                            return hospital.ID,hospital.TravelTime
                        elseif triagecat == 2 && hospital.CapT2 > 0
                            return hospital.ID,hospital.TravelTime
                        elseif triagecat == 3 && hospital.CapT3 > 0
                            return hospital.ID,hospital.TravelTime
                        end
                    else
                        if (triagecat == 1 || triagecat == 4) && hospital.CapT1Ped > 0
                            return hospital.ID,hospital.TravelTime
                        elseif triagecat == 2 && hospital.CapT2Ped > 0
                            return hospital.ID,hospital.TravelTime
                        elseif triagecat == 3 && hospital.CapT3Ped > 0
                            return hospital.ID,hospital.TravelTime
                        end
                    end
                end
            end
        end
        println("No hospital found for $triagecat, $hospitalneed, $adult")
    end
end
function findRoundRobinHospitalId(hospitalQueue ::Array{Hospital,1}, triagecat ::Int64, hospitalneed ::String, adult ::Bool, paramset::Parameterset, inputDB :: SQLite.DB)
    # Returns ID and TravelTime plus puts hospitalID at end of the hospitalOrder for the first hospital in the hospitalOrder list for specified triagecat, facility, and agecategory (adult or peadiatric)
    hospitallist = findHospitalswithFacility(inputDB, hospitalneed, adult)

    for hospitalID in paramset.hospitalOrder
        if length( hospitallist) > 0
            for i in 1:length( hospitallist)
                for hospital in hospitalQueue
                    
                    if hospital.ID == convert(Int32,( hospitallist[i]))
                        if adult
                            if (triagecat == 1 || triagecat == 4) && hospital.CapT1 > 0
                                
                                
                                paramset.hospitalOrder = [deleteat!(paramset.hospitalOrder, [x == hospitalID for x in paramset.hospitalOrder]); hospitalID]
                               
                               
                                return hospital.ID
                            elseif triagecat == 2 && hospital.CapT2 > 0
                               
                               
                                paramset.hospitalOrder = [deleteat!(paramset.hospitalOrder, [x == hospitalID for x in paramset.hospitalOrder]); hospitalID]
                                return hospital.ID
                            elseif triagecat == 3 && hospital.CapT3 > 0
                              
                              
                                paramset.hospitalOrder = [deleteat!(paramset.hospitalOrder, [x == hospitalID for x in paramset.hospitalOrder]); hospitalID]
                                return hospital.ID
                            end
                        else
                            if (triagecat == 1 || triagecat == 4) && hospital.CapT1Ped > 0
                                paramset.hospitalOrder = [deleteat!(paramset.hospitalOrder, [x == hospitalID for x in paramset.hospitalOrder]); hospitalID]
                                return hospital.ID
                            elseif triagecat == 2 && hospital.CapT2Ped > 0
                                paramset.hospitalOrder = [deleteat!(paramset.hospitalOrder, [x == hospitalID for x in paramset.hospitalOrder]); hospitalID]
                                return hospital.ID
                            elseif triagecat == 3 && hospital.CapT3Ped > 0
                                paramset.hospitalOrder = [deleteat!(paramset.hospitalOrder, [x == hospitalID for x in paramset.hospitalOrder]); hospitalID]
                                return hospital.ID
                            end
                        end
                    end
                end
            end
        else
            println("(spreadout) No hospital found for $triagecat, $hospitalneed, $adult")
        end
    end
    
end

function findHospTravelTime(hospitalQueue ::Array{Hospital,1}, hospitalID ::Int,patient::Patient,m,quadrant)
    for hospital in hospitalQueue
        if hospital.ID == hospitalID
                    hospX = hospital.LAT
                    hospY = hospital.LON
                  
                    patientUTM = Geodesy.UTM(patient.x0,patient.y0,0.0)
                    lla_location = LLAfromUTM(quadrant,true,wgs84)
                    patientLLA = lla_location(patientUTM)
                    hospNode = point_to_nodes((hospX,hospY),m)
                    fmpNode = point_to_nodes((patientLLA.lat,patientLLA.lon),m)
                    route, distance,route_time = OpenStreetMapX.fastest_route(m,hospNode,fmpNode)
                  

                    timeToH = (route_time/60.0)*1.0
                   

            return timeToH
            
        end
    end
end
function requestHospital(env ::Environment, hospitalQueue ::Array{Hospital,1}, hospitalID ::Int64, triagecat ::Int64, adult ::Bool)
 
    for hospital in hospitalQueue
      
        if hospital.ID == hospitalID
         
            if adult
                if (triagecat == 1 || triagecat == 4)
                    if hospital.CapT1 > 0
                    hospital.CapT1 -= 1
                        
                    else
                        return
                    end
              
                    if hospital.CapT1 <= 0
                        
                        @process checkForSurge(env, hospital)
                    end
                    return
                elseif triagecat == 2
                    hospital.CapT2 -= 1
           
                    if hospital.CapT2 <= 0
                        @process checkForSurge(env, hospital)
                    end
                    return
                elseif triagecat == 3
                    hospital.CapT3 -= 1
             
                    if hospital.CapT3 <= 0
                        @process checkForSurge(env, hospital)
                    end
                    return
                end
            else
                if (triagecat == 1 || triagecat == 4)
                    hospital.CapT1Ped -= 1
                    if hospital.CapT1Ped <= 0
                        @process checkForSurge(env, hospital)
                    end
                    return
                elseif triagecat == 2
                    hospital.CapT2Ped -= 1
                    if hospital.CapT2Ped <= 0
                        @process checkForSurge(env, hospital)
                    end
                    return
                elseif triagecat == 3
                    hospital.CapT3Ped -= 1
                    if hospital.CapT3Ped <= 0
                        @process checkForSurge(env, hospital)
                    end
                    return
                end
            end
        
        end
    end
    
    error("Hospital $hospitalID not found in the hospitalqueue")
end
function checkSupervision(triagecat ::Int64, supervisionlevels ::Array, supervisionavailable ::Tuple{Bool,Bool})
    
    if triagecat == 9
        return false
    end
    if triagecat == 4
        triagecat = 1
    end
    if triagecat == 5
        triagecat = 3
    end
    if triagecat == 0
        triagecat = 1
    end
    supmedlevel = supervisionlevels[triagecat] # allowed supervision
    
    mediclevel = 0

    
    for j in 1:lastindex(supmedlevel)
        mediclevel = supmedlevel[j]
       
        if (mediclevel == 0  && supervisionavailable[1] && supervisionavailable[2]) 
            return mediclevel
        end
        if (mediclevel == 1 && supervisionavailable[1])
            return mediclevel
        end
        if (mediclevel == 2 )
           
            return mediclevel
        end
        if mediclevel == 3
            return mediclevel
        end
    end

    return false

end
function findVicPosition(score ::Int64, queue :: Array)
    if length(queue) >= 0
        for i = 1:length(queue)
            qscore = queue[i][1]
            if qscore > score
                index = i
                return index
            end
        end
    end
        index = length(queue)+1
   
        return index
    
end
function enterVicTranspQueue!(queue :: Array, patient ::Patient, score ::Int64)
    # put victim in queue in correct position (based on score)
    
    index = findVicPosition(score, queue)

    if index <= length(queue)
        
           
      insert!(queue,index,(score,patient))
     
    else
      
       push!(queue,(score,patient))
          
    end
end
@resumable function checkForSurge(env ::Environment, hospital ::Hospital)
    if !hospital.Surgemode
        hospital.Surgemode = true
      
        @process updateCapacity(env, hospital)
        
    end
end
@resumable function updateCapacity(env ::Environment, hospital ::Hospital)
    if hospital.CapT1 < hospital.SurgeT1
        hospital.CapT1 = hospital.SurgeT1
    end
    if hospital.CapT2 < hospital.SurgeT2
        hospital.CapT2 = hospital.SurgeT2
    end
    if hospital.CapT3 < hospital.SurgeT3
        hospital.CapT3 = hospital.SurgeT3
    end

    if hospital.CapT1Ped < hospital.SurgeT1Ped
        hospital.CapT1Ped = hospital.SurgeT1Ped
    end
    if hospital.CapT2Ped < hospital.SurgeT2Ped
        hospital.CapT2Ped = hospital.SurgeT2Ped
    end
    if hospital.CapT3Ped < hospital.SurgeT3Ped
        hospital.CapT3Ped = hospital.SurgeT3Ped
    end

    @yield timeout(env, 60.0)
    
  

end
function stopSurge(hospitalQueue ::Array{Hospital,1})
    for hospital in hospitalQueue
        hospital.Surgemode = false
    end
end
function checkMobility(env::Environment,patient::Patient)
    #convert mobility int to bool for if purposes
    canstand = patient.mobile
    
    return (canstand == 1 && patient.incapacitated == 0)
     
end
function calculateISS(firstInjury::Injury,secondInjury::Injury,thirdInjury::Injury)
    location1 = firstInjury.injuryLocation
    location2 = secondInjury.injuryLocation
    location3 = thirdInjury.injuryLocation

    if (location1 != location2) && (location1 != location3) && (location2 != location3)
        iss = firstInjury.severityScore^2 + secondInjury.severityScore^2 + thirdInjury.severityScore^2
    elseif (location1 == location2) && (location1 != location3) 
        iss = max(firstInjury.severityScore^2 ,secondInjury.severityScore^2) + thirdInjury.severityScore^2
    elseif  (location1 == location3) && (location2 != location1) 
        iss = max(firstInjury.severityScore^2 ,thirdInjury.severityScore^2) + secondInjury.severityScore^2 
    else
        iss = max(firstInjury.severityScore^2 ,thirdInjury.severityScore^2,secondInjury.severityScore^2)
        
    end
    return iss

end
function generateVictimPosition(env::Environment,  event::Event,location::UTM,consts::Constants)
   
    
    
    
    positions =zeros(Float64,event.casualtyEstimate,8)
    k=0 #from 0 to numberofVictims
    if event.sustained == true
       
        spread = event.aoe*5 #assuming 200x200m spread (look up)
        for j in 1:trunc(Int32,event.threatLength)
            hitZoneX = spread*(-1.0+rand()*2)
            hitZoneY = spread*(-1.0+rand()*2)
        
            
            for i in 1:event.casualtyRate #(victims per hit)
                
                randDistX= (0.5 + rand()*(event.aoe)*2*sin(2*3.14159*rand())) #spawn victims at a x of 0.5+-*AOE*4m around the event epicenter
                randDistY= (0.5 + rand()*(event.aoe)*2*sin(2*3.14159*rand())) #spawn victims at a y of 0.5+-*AOE*4m around the event epicenter
               
                k=k+1
               
                positions[k,1] =location.x+hitZoneX+randDistX #victim position X UTM
                positions[k,2] =location.y+hitZoneY+randDistY #victim position Y UTM
                vic_utm = Geodesy.UTM(positions[k,1],positions[k,2],0.0) 
                lla_location = LLAfromUTM(consts.quadrant,true,wgs84)
                vic_lla=lla_location(vic_utm)
             
                positions[k,3] = location.x+hitZoneX
                positions[k,4] = location.y+hitZoneY
                hit_utm = Geodesy.UTM(positions[k,3],positions[k,4],0.0)
                hit_lla=lla_location(hit_utm)
                positions[k,5] = vic_lla.lat
                positions[k,6] = vic_lla.lon
                positions[k,7] = hit_lla.lat
                positions[k,8] = hit_lla.lon
                               
            end
        
        end
        
    else
        for i in 1:event.casualtyEstimate
            hitZoneX = 0.0
            hitZoneY = 0.0
            randDistX= (0.5 + rand()*(event.aoe)*2*sin(2*3.14159*rand())) #spawn victims at a x of 0-*AOE*m around the event epicenter
            randDistY= (0.5 + rand()*(event.aoe)*2*sin(2*3.14159*rand())) #spawn victims at a y of 0-*AOE*m around the event epicenter
            positions[i,1] =location.x+hitZoneX+randDistX #victim position X
            positions[i,2] =location.y+hitZoneY+randDistY #victim position Y
            vic_utm = Geodesy.UTM(positions[i,1],positions[i,2],0.0) 
            lla_location = LLAfromUTM(consts.quadrant,true,wgs84)
            vic_lla=lla_location(vic_utm)
            
            positions[i,3] = location.x+hitZoneX
            positions[i,4] = location.y+hitZoneY
            hit_utm = Geodesy.UTM(positions[i,3],positions[i,4],0.0)
            hit_lla=lla_location(hit_utm)
            positions[i,5] = vic_lla.lat
            positions[i,6] = vic_lla.lon
            positions[i,7] = hit_lla.lat
            positions[i,8] = hit_lla.lon
        end
    end
    return positions
end
function generateVictimPositionRad(env::Environment,  event::Event,location::UTM,consts::Constants)
   
    
    
    positions =zeros(Float64,consts.numberofVictims,8)
   
    positions =zeros(Float64,consts.numberofVictims,8)
   
    k=0 #from 0 to numberofVictims
    if event.sustained == true
        
        for j in 1:trunc(Int32,event.threatLength)
            hitZoneX = 0.0
            hitZoneY = 0.0
        
            
            for i in 1:event.casualtyRate #(victims per hit)
                
              
                randDistX= (0.5 + rand()*(event.aoe)*25000*sin(2*3.14159*rand())) #spawn victims at a x of 0.5+-*AOE*4m around the event epicenter
                randDistY= (0.5 + rand()*(event.aoe)*25000*sin(2*3.14159*rand())) #spawn victims at a y of 0.5+-*AOE*4m around the event epicenter
                k=k+1
                positions[k,1] =location.x+hitZoneX+randDistX #victim position X
                positions[k,2] =location.y+hitZoneY+randDistY #victim position Y
                vic_utm = Geodesy.UTM(positions[k,1],positions[k,2],0.0) 
                lla_location = LLAfromUTM(consts.quadrant,true,wgs84)
                vic_lla=lla_location(vic_utm)
             
                positions[k,3] = location.x+hitZoneX
                positions[k,4] = location.y+hitZoneY
                hit_utm = Geodesy.UTM(positions[k,3],positions[k,4],0.0)
                hit_lla=lla_location(hit_utm)
                positions[k,5] = vic_lla.lat
                positions[k,6] = vic_lla.lon
                positions[k,7] = hit_lla.lat
                positions[k,8] = hit_lla.lon
                               
            end
        
        end
        
    else #tacNuke
        for i in 1:consts.numberofVictims
            hitZoneX = 0.0
            hitZoneY = 0.0
            x,y = rand(1600.0:5000.0)*randn(2)
            while sqrt(x^2 + y^2) < 1600.0 || sqrt(x^2 + y^2) > 5000.0
                x,y = rand(1600.0:5000.0)*randn(2)
            end
            positions[i,1], positions[i,2] = x + location.x,y + location.y
            # randDistX= min(10000,(10000 + rand()*(event.aoe)*sin(1/2*3.14159*rand())))*rand([-1:1]) #spawn victims past a circle of 
            # randDistY= sqrt(10000^2-randDistX^2)
            # positions[i,1] =location.x+hitZoneX+randDistX #victim position X
            # positions[i,2] =location.y+hitZoneY+randDistY #victim position Y
            vic_utm = Geodesy.UTM(positions[i,1],positions[i,2],0.0) 
            lla_location = LLAfromUTM(consts.quadrant,true,wgs84)
            vic_lla=lla_location(vic_utm)
            
            positions[i,3] = location.x+hitZoneX
            positions[i,4] = location.y+hitZoneY
            hit_utm = Geodesy.UTM(positions[i,3],positions[i,4],0.0)
            hit_lla=lla_location(hit_utm)
            positions[i,5] = vic_lla.lat
            positions[i,6] = vic_lla.lon
            positions[i,7] = hit_lla.lat
            positions[i,8] = hit_lla.lon
        end
    end
    return positions
end
function walkable(r::Float64,area::Float64) #incapacity based on Dullum2010 area is the lethal area

    prob = exp(-pi*r^2/area)
    t= rand()
    if t < prob
        walk = 0
    else 
        walk = 1
    end
    return walk,t,prob
end
function shrapnelPenetration(r::Float64,nFragments::Int64)

    return 1-exp(-nFragments*0.5/(4*pi*r^2)) #dullum2010
end
function updateHospitalNeed(patient::Patient)
    #Find hospital with that fits all requirements
  
    maxISS = patient.iss
    
    if maxISS == 1 || maxISS == 0
        if patient.hospitalneed == "H"
            patient.hospitalneed = "H"
        else
            patient.hospitalneed = patient.hospitalneed
        end
    end
    if maxISS == 2
       if patient.hospitalneed == "H"
            patient.hospitalneed = "OPC"
       else
            patient.hospitalneed =patient.hospitalneed
       end
    end
    if patient.iss > 25
        patient.hospitalneed = "PT"
    end
    if maxISS >= 3
        patient.hospitalneed = "S"
        if patient.victimInjury.injuryLocation == "head/neck" || patient.victimInjury.injuryLocation == "spine"
            patient.hospitalneed = "NS"
        
        end
        if patient.victimInjury.injuryLocation == "chest"
            patient.hospitalneed = "CTS"
           
        end
    end
    patient.hospitalneed = patient.hospitalneed
end
function scoutForMedic(pos:: Tuple{Float64,Float64},radius::Float64)

    #returns a list of victimIDS who are within radius of the victim with pos(x,y) and are able to provide help
    
    Patientlist =DataFrame(CSV.File("victimsAll.txt"))
    for i in 1:size(Patientlist,1)
        distTopatient = sqrt((Patientlist[i,7]-pos[1])^2 +(Patientlist[i,8]-pos[2])^2)
        if distTopatient < 25.0 &&  Patientlist[i,15]==1
      
            return 1,i
        end
    end
    return 0,0
    

end
function getDistance(start::UTM,destination::UTM)
    
    dist = sqrt((start.x - destination.x)^2 + (start.y - destination.y)^2)
   
    return dist
end
function randomSort()

    randomsort = rand(1:100)
   
    return randomsort

end
function findVictimTransport(env::Environment,transportqueue ::Array, hospitalQueue ::Array{Hospital,1}, supervisionlevels ::Array, supervisionavailable ::Tuple{Bool,Bool}, HospitalDistribution ::String,inputDB,paramset::Parameterset,locationx,locationy,m,quadrant,offroad)
    if length(transportqueue) == 0
      
        return nothing,nothing,nothing
    end
    i = 1
    continuecheck = true
      
    while continuecheck
        if i > length(transportqueue)
            return nothing,nothing,nothing
        end
        victim
        (score, patient) = transportqueue[i]
        
        if isdead(patient)
            
            deleteat!(transportqueue,i)
            return
        end
        if paramset.Policy == "StayPlay" 
            # correct triage after treatment
           if patient.triage !=9
                triagecat = patient.triage
               
           else triagecat =9
                return nothing,nothing,nothing
           end
        elseif paramset.Policy == "ScoopRun" 
            # incorrect triage
            if patient.triage !=9
                triagecat = triageMistakes(patient.triage, patient.triagecorrectness)
           else triagecat =9
                return nothing,nothing,nothing
           end
                                
        end
      
        locationx=patient.x0
        locationy=patient.y0
        supervision = checkSupervision(triagecat, supervisionlevels, supervisionavailable)
            
        if supervision === false # no supervision available, need 3*= otherwise 0 == false -> true
            
            i += 1 # check next victim
            # else, supervision available, check hospital
      
        elseif checkTotalCapacity(hospitalQueue, triagecat, patient.hospitalneed, patient.age >= 15,inputDB) <= 0 # no hospital available
           
            i += 1 # check next victim
            
        else
            # all checks positive, remove victim from queue, return victim with supervisionlevel
            if HospitalDistribution == "CloseFirst"
              
                hospid,t = findClosestHospitalId(hospitalQueue, triagecat, patient.hospitalneed, patient.age >= 15,inputDB,locationx,locationy,m,quadrant,offroad)
                patient.hospitalID = hospid
                
                
                            
            elseif HospitalDistribution == "SpreadOut"
                hospid = findRoundRobinHospitalId(hospitalQueue, triagecat, patient.hospitalneed, patient.age >= 15, paramset,inputDB)
                patient.hospitalID = hospid
                
                
            end
        
            if length(transportqueue) > 0
           
                deleteat!(transportqueue,i)
            end
         
            return patient, supervision, hospid
        end
       
    end
end
@resumable function createThreat(env::Environment, threatNumber ::Int64,scheduledTime::Float64,epicenter ::UTM,duration::Float64,outputDB ::SQLite.DB, firefighters ::Resource, boatList::BoatList,mediclist ::MedicList, ambulist ::AmbuList,medevac::Medevac, transportqueue ::Array, transportqueue1::Array,transportqueueT3 ::Array, hospitalQueue ::Array, paramset ::Parameterset,SS_array::Array,treatment_data::Array,consts::Constants,inputDB,m,fmp_LLA,ccp_LLA,homeUTM,ccpUTM,timeCCPHome,timeccpFMP,timeFMPHome,fmpUTM,threatId,offroad)

    if threatNumber == 1
        name = "artillery155mm"
        roe = 40.0 # radius of effect up to 500 meters following champion2009
        schedule = scheduledTime
        magnitude = 3
        location = epicenter
        cbrne = 5 # it is an explosion
        threatLength = duration #option Artillery strike for *duration* minutes
        casualtyRate = 13 #13 victims/hit
        consts.numberofVictims = trunc(Int64,duration*casualtyRate)
        threat = Event(name,roe,schedule,magnitude,epicenter,cbrne,threatLength,true,casualtyRate,consts)
        
    end
    if threatNumber == 2
        name = "tacNuke"
        aoe = 1.0
        schedule = scheduledTime
        magnitude = 3
        location = epicenter
        cbrne = 4 # it is a nuclear threat
        threatLength = duration #instantaneous
        casualtyRate = 1000 #NA because it is instant
        consts.numberofVictims = trunc(Int64,duration*casualtyRate) 
        threat = Event(name,aoe,schedule,magnitude,epicenter,cbrne,threatLength,false,casualtyRate,consts)
        
    end
    if threatNumber == 3
        name = "Chem" #demo
        aoe = 4.0
        schedule = scheduledTime
        magnitude = 2
        location = epicenter
        cbrne = 5 # it is an explosion
        consts.numberofVictims = trunc(Int64,duration*casualtyRate) #person in front
        casualtyRate = 50 #NA because it is instant
        consts.numberofVictims = trunc(Int64,duration*casualtyRate) #person in front
        threat = Event(name,aoe,schedule,magnitude,epicenter,cbrne,threatLength,false,casualtyRate,consts)
    end
    if threatNumber == 4
        name = "sustained762_fire"
        aoe = 1.0
        schedule = scheduledTime #no threat for armored vehicles
        magnitude = 1
        location = epicenter
        cbrne = 5 # it is an explosion
        threatLength = duration #continous or instant
        casualtyRate = 3 #NA function of the armor and placement!
        consts.numberofVictims = trunc(Int64,duration*casualtyRate) #person in front
        threat = Event(name,aoe,schedule,magnitude,epicenter,cbrne,threatLength,true,casualtyRate,consts)
    end
    if threatNumber == 5
        name = "FPVDrone"
        roe = 5.0
        schedule = scheduledTime
        magnitude = 1
        location = epicenter
        cbrne = 5 # it is an explosion
        threatLength = duration #instantaneous
        casualtyRate = 8 #NA because it is instant
        consts.numberofVictims = trunc(Int64,duration*casualtyRate) #person in front
        threat = Event(name,roe,schedule,magnitude,epicenter,cbrne,threatLength,false,casualtyRate,consts)
    end
    if threatNumber == 6
        name = "droneStrike"
        roe = 10.0 # radius of effect 
        schedule = scheduledTime
        magnitude = 3
        location = epicenter
        cbrne = 5 # it is an explosion
        threatLength = duration #option Artillery strike for *duration* minutes
        casualtyRate = 30 #30 victims/hit
        consts.numberofVictims = trunc(Int64,duration*casualtyRate)
        threat = Event(name,roe,schedule,magnitude,epicenter,cbrne,threatLength,false,casualtyRate,consts)
        
    end
    if threatNumber == 7
        name = "RefugeeEvacuation"
        roe = 20.0 # radius of effect 
        schedule = scheduledTime
        magnitude = 1
        location = epicenter
        cbrne = 0 # it is an explosion
        threatLength = duration #option 
        casualtyRate = 400 #5 victims/cluster
        consts.numberofVictims = trunc(Int64,duration*casualtyRate)
        threat = Event(name,roe,schedule,magnitude,epicenter,cbrne,threatLength,false,casualtyRate,consts)
        
    end
    if threatNumber == 8
        name = "TIC"
        roe = 10.0 # radius of effect 
        schedule = scheduledTime
        magnitude = 3
        location = epicenter
        cbrne = 5 # it is an explosion
        threatLength = duration #option Artillery strike for *duration* minutes
        casualtyRate = 4 #4 victims/hit
        consts.numberofVictims = trunc(Int64,duration*casualtyRate)
        threat = Event(name,roe,schedule,magnitude,epicenter,cbrne,threatLength,false,casualtyRate,consts)
        
    end
    if threatNumber == 6
        name = "droneStrike"
        roe = 10.0 # radius of effect 
        schedule = scheduledTime
        magnitude = 3
        location = epicenter
        cbrne = 5 # it is an explosion
        threatLength = duration #option Artillery strike for *duration* minutes
        casualtyRate = 30 #30 victims/hit
        consts.numberofVictims = trunc(Int64,duration*casualtyRate)
        threat = Event(name,roe,schedule,magnitude,epicenter,cbrne,threatLength,false,casualtyRate,consts)
        
    end
    if threatNumber == 7
        name = "RefugeeEvacuation"
        roe = 20.0 # radius of effect 
        schedule = scheduledTime
        magnitude = 1
        location = epicenter
        cbrne = 0 # it is an explosion
        threatLength = duration #option 
        casualtyRate = 400 #5 victims/cluster
        consts.numberofVictims = trunc(Int64,duration*casualtyRate)
        threat = Event(name,roe,schedule,magnitude,epicenter,cbrne,threatLength,false,casualtyRate,consts)
        
    end
    if threatNumber == 8
        name = "TIC"
        roe = 10.0 # radius of effect 
        schedule = scheduledTime
        magnitude = 3
        location = epicenter
        cbrne = 5 # it is an explosion
        threatLength = duration #option Artillery strike for *duration* minutes
        casualtyRate = 4 #4 victims/hit
        consts.numberofVictims = trunc(Int64,duration*casualtyRate)
        threat = Event(name,roe,schedule,magnitude,epicenter,cbrne,threatLength,false,casualtyRate,consts)
        
    end

       @process spawnVictims(env, threat,SS_array, treatment_data,threat.location,outputDB, firefighters, boatList,mediclist, ambulist,medevac, transportqueue,transportqueue1, transportqueueT3, hospitalQueue, paramset,consts,inputDB,m,fmp_LLA,ccp_LLA,homeUTM,ccpUTM,timeCCPHome,timeccpFMP,timeFMPHome,fmpUTM,threatId,offroad)
end
@resumable function readVictims(env ::Environment, numberofVictims::Int64, firefighters ::Resource, mediclist ::MedicList, ambulist ::AmbuList,medevac::Medevac, transportqueue ::Array, transportqueueT3 ::Array, hospitalQueue ::Array, paramset ::Parameterset,SS_array::Array,treatment_data::Array,consts::Constants,outputDB,homeUTM,ccpUTM,fmp_LLA,ccp_LLA,timeCCPHome,timeccpFMP,timeFMPHome,fmpUTM,inputDB,m,cbrn,disaster)
   
    unit = DeconUnit(env,20,1000,40,100,1000,300,1000,100)
    amst = Amst(env,100,40,false)
    if cbrn == 5
        Patientslist =DataFrame(CSV.File("victimsArtillery.txt"))
       # Patientslist =DataFrame(CSV.File("victimsAll.txt"))
    end
    if cbrn == 4
        Patientslist =DataFrame(CSV.File("victimsRad.txt"))
    end
   nrow(Patientslist)
    
    for i in 1:consts.numberofVictims #generate numberofVictims victims
        if cbrn == 4
            victag = Patientslist[i,1]#shift caused by chemical victims 
        
            age = Patientslist[i,2]
            triage = Patientslist[i,3]
            is_treated = 0
            canwalk = Patientslist[i,15]
            
            x0 = Patientslist[i,7]
            y0 = Patientslist[i,8]
            
            
            iss = Patientslist[i,10]
            timeOfDeath = min(Patientslist[i,14],Patientslist[i,20])

            hitloc = (Patientslist[i,5],Patientslist[i,6])
            distFromBlast = Patientslist[i,9]
            armor = 100.0
            bleed = Patientslist[i,16]
        
            timeOfDeath = convert(Float64,timeOfDeath)
            ip = 7
        
            SimParams = [Patientslist[i,11]  Patientslist[i,12]  Patientslist[i,13]  timeOfDeath ip Patientslist[i,17]  Patientslist[i,18]  Patientslist[i,19] Patientslist[i,20]]
            if (Patientslist[i,13] == 1.0 || Patientslist[i,19] == 1)
            lethally_injured = 1
            else
            lethally_injured = 0
            end
            maxSimScore = 20.0
            treatmentParams=[0.0,0.0,0.0,0.0,0.0]
            
            injury =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            injury2 =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            injury3 =Injury(1,"","","",1,[""],"",6.0,"",0.0)

            injuryn = 0
            injuryn2= 0
            injuryn3 = 0
            
            injury = milInjurySetter(injuryn,injury)
            injury2 = milInjurySetter(injuryn2,injury2)
            injury3 = milInjurySetter(injuryn3,injury3)
            
                hospitalNeed = injury.hospitalNeed
                if triage == 5
                    hospitalNeed = "H"
                    injury.hospitalNeed = "H"
                    injury2.hospitalNeed = "H"
                    injury3.hospitalNeed = "H"
                end

            if paramset.TriageLevel == "None"
                triage = 0
            end
            patient = Patient(env, victag, age, is_treated, canwalk , lethally_injured,triage,x0,y0,iss,timeOfDeath,0.0,SimParams,maxSimScore,bleed,(0,0),injury,injury2,injury3,hospitalNeed,consts,SS_array,distFromBlast,treatmentParams,0,"PoI",ip)
            
                
            patient.hospitalneed=updateHospitalNeed(patient)
            sscore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
                    
            if sscore <= 0.01 && patient.isdead == 0

                patient.triage = 9
                patient.isdead = 1
                VictimLogOutputDB(outputDB, patient.vicid,  now(env), 188,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,0,patient.SimParams[3],patient.bloodvolume)
                
            end
            if patient.triage !=9
            if consts.logg
                VictimLogOutputDB(outputDB, patient.vicid, now(env), 42,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,0,patient.SimParams[3],patient.bloodvolume)
            end
            end
        # push!(victims,patient)
        #  for t in 1:1000
                
        #        time_from_onset = convert(Float64,t*60)
        #        SimScore = generateSimedisScoreGompertz(env,time_from_onset,patient.SimParams)
        #        displayTime = time_from_onset/60.
        #        SS_array = push!(SS_array,[patient.vicid  displayTime  SimScore])
        #  end

        @yield request(patient.transport) #resource to let the victim function continue after transport
        @process victim(env, patient, outputDB, firefighters, unit,amst,mediclist, ambulist,medevac, transportqueue, transportqueueT3, hospitalQueue, paramset,treatment_data,fmpUTM,homeUTM,ccpUTM,fmp_LLA,ccp_LLA,timeCCPHome,timeccpFMP,timeFMPHome,consts,inputDB,m,offroad)
                
        end
        if cbrn==5
            victag = Patientslist[i,1]#shift caused by chemical victims 
        
            age = Patientslist[i,2]
            triage = Patientslist[i,3]
            is_treated = 0
            canwalk = Patientslist[i,15]
            
            x0 = Patientslist[i,7]
            y0 = Patientslist[i,8]
            
            iss = Patientslist[i,10]
            timeOfDeath = Patientslist[i,14]

            hitloc = (Patientslist[i,5],Patientslist[i,6])
            distFromBlast = Patientslist[i,9]
            armor = 100.0
            bleed = Patientslist[i,16]
        
            hitloc = (Patientslist[i,5],Patientslist[i,6])
            distFromBlast = Patientslist[i,9]
            armor = 100.0
            bleed = Patientslist[i,16]
        

            
            timeOfDeath = convert(Float64,timeOfDeath)
            ip = 7
        
            SimParams = [Patientslist[i,11]  Patientslist[i,12]  Patientslist[i,13]  timeOfDeath ip 0.0  0.0  0.0 0.0]
            if (Patientslist[i,13] == 1.0)
            lethally_injured = 1
            else
            lethally_injured = 0
            end
            maxSimScore = 20.0
            treatmentParams=[0.0,0.0,0.0,0.0,0.0]
            
            injury =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            injury2 =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            injury3 =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            
            timeOfDeath = convert(Float64,timeOfDeath)
            ip = 7
        
            SimParams = [Patientslist[i,11]  Patientslist[i,12]  Patientslist[i,13]  timeOfDeath ip 0.0  0.0  0.0 0.0]
            if (Patientslist[i,13] == 1.0)
            lethally_injured = 1
            else
            lethally_injured = 0
            end
            maxSimScore = 20.0
            treatmentParams=[0.0,0.0,0.0,0.0,0.0]
            
            injury =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            injury2 =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            injury3 =Injury(1,"","","",1,[""],"",6.0,"",0.0)

            injuryn = 0
            injuryn2= 0
            injuryn3 = 0
            
            injury = milInjurySetter(injuryn,injury)
            injury2 = milInjurySetter(injuryn2,injury2)
            injury3 = milInjurySetter(injuryn3,injury3)
            
                hospitalNeed = injury.hospitalNeed
                if triage == 5
                    hospitalNeed = "H"
                    injury.hospitalNeed = "H"
                    injury2.hospitalNeed = "H"
                    injury3.hospitalNeed = "H"
                end
            injuryn = 0
            injuryn2= 0
            injuryn3 = 0
            
            injury = milInjurySetter(injuryn,injury)
            injury2 = milInjurySetter(injuryn2,injury2)
            injury3 = milInjurySetter(injuryn3,injury3)
            
                hospitalNeed = injury.hospitalNeed
                if triage == 5
                    hospitalNeed = "H"
                    injury.hospitalNeed = "H"
                    injury2.hospitalNeed = "H"
                    injury3.hospitalNeed = "H"
                end

            if paramset.TriageLevel == "None"
                triage = 0
            end
            
            
                patient = Patient(env, victag, age, is_treated, canwalk , lethally_injured,triage,x0,y0,iss,timeOfDeath,0.0,SimParams,maxSimScore,bleed,(0,0),injury,injury2,injury3,hospitalNeed,consts,SS_array,distFromBlast,treatmentParams,0,"PoI",ip)
            
                
                patient.hospitalneed=updateHospitalNeed(patient)
                sscore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
                        
                if sscore <= 0.01 && patient.isdead == 0

                    patient.triage = 9
                    patient.isdead = 1
                    VictimLogOutputDB(outputDB, patient.vicid,  now(env), 188,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,0,patient.SimParams[3],patient.bloodvolume)
                    
                end
                if patient.triage !=9
                if consts.logg
                    VictimLogOutputDB(outputDB, patient.vicid, now(env), 42,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,0,patient.SimParams[3],patient.bloodvolume)
                end
                
        

                @yield request(patient.transport) #resource to let the victim function continue after transport
                @process victim(env, patient, outputDB, firefighters,unit,amst, mediclist, ambulist,medevac, transportqueue, transportqueueT3, hospitalQueue, paramset,treatment_data,fmpUTM,homeUTM,ccpUTM,fmp_LLA,ccp_LLA,timeCCPHome,timeccpFMP,timeFMPHome,consts,inputDB,m,offroad)
        
            end
          
        end
    end

    
    @process lastVictim(env, firefighters, mediclist, ambulist,medevac, paramset.Policy,paramset,outputDB,consts)
  
   
end
@resumable function readVictimsMulti(env ::Environment, firefighters ::Resource, mediclist ::MedicList, ambulist ::AmbuList,medevac::Medevac, transportqueue ::Array, transportqueueT3 ::Array, hospitalQueue ::Array, paramset ::Parameterset,SS_array::Array,treatment_data::Array,consts::Constants,outputDB,fmpUTM,homeUTM,ccpUTM,fmp_LLA,ccp_LLA,timeCCPHome,timeccpFMP,timeFMPHome,inputDB,m,cbrn,disaster,offroad)
   
    unit = DeconUnit(env,20,1000,40,100,1000,300,1000,100)
    amst = Amst(env,100,40,false)
    if cbrn == 5

        Patientslist =DataFrame(CSV.File("victimsAll.txt"))
    end
    if cbrn == 4
        Patientslist =DataFrame(CSV.File("victimsRad.txt"))
    end
    #@show nrow(Patientslist)

    for i in 1:size(Patientslist,1) #generate numberofVictims victims (victims per threat times number of threats)
    
        if cbrn == 4
            victag = Patientslist[i,1]#shift caused by chemical victims 
        
            age = Patientslist[i,2]
            triage = Patientslist[i,3]
            is_treated = 0
            canwalk = Patientslist[i,15]
            
           
            patient_LLA = Geodesy.LLA(Patientslist[i,17],Patientslist[i,18]) #LLA coord for patient
            utm_location = UTMfromLLA(consts.quadrant,true,wgs84)
            patientUTM = utm_location(patient_LLA)
            x0 = patientUTM.x
            y0= patientUTM.y
            
            iss = Patientslist[i,10]
            timeOfDeath = min(Patientslist[i,14],Patientslist[i,20])

            hitloc = (Patientslist[i,5],Patientslist[i,6])
            distFromBlast = Patientslist[i,9]
            bloodvolume = 5.0
            bleed = Patientslist[i,16]
        
            timeOfDeath = convert(Float64,timeOfDeath)
            ip = 7
        
            SimParams = [Patientslist[i,11]  Patientslist[i,12]  Patientslist[i,13]  timeOfDeath ip Patientslist[i,17]  Patientslist[i,18]  Patientslist[i,19] Patientslist[i,20]]
            if (Patientslist[i,13] == 1.0 || Patientslist[i,19] == 1)
                lethally_injured = 1
            else
                lethally_injured = 0
            end
            maxSimScore = 20.0
            treatmentParams=[0.0,0.0,0.0,0.0,0.0]
            
            injury =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            injury2 =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            injury3 =Injury(1,"","","",1,[""],"",6.0,"",0.0)

            injuryn = 0
            injuryn2= 0
            injuryn3 = 0
            
            injury = milInjurySetter(injuryn,injury)
            injury2 = milInjurySetter(injuryn2,injury2)
            injury3 = milInjurySetter(injuryn3,injury3)
            
                hospitalNeed = injury.hospitalNeed
                if triage == 5
                    hospitalNeed = "H"
                    injury.hospitalNeed = "H"
                    injury2.hospitalNeed = "H"
                    injury3.hospitalNeed = "H"
                end

            if paramset.TriageLevel == "None"
                triage = 0
            end
            scheduledTime = 0.0
            facility = convert(String,Patientslist[i,23])
            facilities = ["PEACH","APPLE","PEAR","WALNUT","PECAN","RAISIN","BANANA","BLACKBERRY","ORANGE","PLUM","MANGO","PINEAPPLE","CHESTNUT"]

            if facility in facilities
                facility = "PoI"
            end
            patient = Patient(env, victag, age, is_treated, canwalk , lethally_injured,triage,x0,y0,iss,timeOfDeath,scheduledTime,SimParams,maxSimScore,bleed,(0,0),injury,injury2,injury3,hospitalNeed,consts,SS_array,distFromBlast,treatmentParams,0,facility,ip)
            
                
            patient.hospitalneed=updateHospitalNeed(patient)
            @yield timeout(env,scheduledTime)
            sscore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
                    
            if sscore <= 0.01 && patient.isdead == 0

                patient.triage = 9
                patient.isdead = 1
                VictimLogOutputDB(outputDB, patient.vicid,  now(env), 188,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,0,patient.SimParams[3],patient.bloodvolume)
                
            end
            if patient.triage !=9
                if consts.logg
                    VictimLogOutputDB(outputDB, patient.vicid, now(env), 42,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,0,patient.SimParams[3],patient.bloodvolume)
                end
            end
            # push!(victims,patient)
            #  for t in 1:1000
                    
            #        time_from_onset = convert(Float64,t*60)
            #        SimScore = generateSimedisScoreGompertz(env,time_from_onset,patient.SimParams)
            #        displayTime = time_from_onset/60.
            #        SS_array = push!(SS_array,[patient.vicid  displayTime  SimScore])
            #  end
              
            if patient.triage !=9
                if consts.logg
                    VictimLogOutputDB(outputDB, patient.vicid, now(env), 42,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,0,patient.SimParams[3],patient.bloodvolume)
                end
            end
            # push!(victims,patient)
            #  for t in 1:1000
                    
            #        time_from_onset = convert(Float64,t*60)
            #        SimScore = generateSimedisScoreGompertz(env,time_from_onset,patient.SimParams)
            #        displayTime = time_from_onset/60.
            #        SS_array = push!(SS_array,[patient.vicid  displayTime  SimScore])
            #  end

            @yield request(patient.transport) #resource to let the victim function continue after transport
            @process victim(env, patient, outputDB, firefighters,unit,amst, mediclist, ambulist,medevac, transportqueue, transportqueueT3, hospitalQueue, paramset,treatment_data,fmpUTM,homeUTM,ccpUTM,fmp_LLA,ccp_LLA,timeCCPHome,timeccpFMP,timeFMPHome,consts,inputDB,m,offroad)
                
        end
        
        if cbrn==5
            victag = Patientslist[i,1]#
            
        
            age = Patientslist[i,2]
            triage = Patientslist[i,3]
            is_treated = 0
            canwalk = Patientslist[i,15]
            
            patient_LLA = Geodesy.LLA(Patientslist[i,17],Patientslist[i,18]) #LLA coord for patient
            utm_location = UTMfromLLA(consts.quadrant,true,wgs84)
            patientUTM = utm_location(patient_LLA)
          
            x0 = patientUTM.x
        
            y0= patientUTM.y
            
            
            iss = Patientslist[i,10]
            timeOfDeath = Patientslist[i,14]

            hitloc = (Patientslist[i,5],Patientslist[i,6])
            distFromBlast = Patientslist[i,9]
            bloodvolume = 5.0
            bleed = Patientslist[i,16]
        

            
            timeOfDeath = convert(Float64,timeOfDeath)
            scheduledTime = convert(Float64,Patientslist[i,21])
            if consts.inputDBname=="vw24.sqlite"
                ip = Patientslist[i,25]
            else 
                ip = 7
            end
        
            SimParams = [Patientslist[i,11]  Patientslist[i,12]  Patientslist[i,13]  timeOfDeath ip 0.0  0.0  0.0 0.0]
           
            if (Patientslist[i,13] == 1.0)
            lethally_injured = 1
            else
            lethally_injured = 0
            end
            maxSimScore = 20.0
            treatmentParams=[0.0,0.0,0.0,0.0,0.0]
            
            injury =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            injury2 =Injury(1,"","","",1,[""],"",6.0,"",0.0)
            injury3 =Injury(1,"","","",1,[""],"",6.0,"",0.0)

            injuryn = 0
            injuryn2= 0
            injuryn3 = 0
            
            injury = milInjurySetter(injuryn,injury)
            injury2 = milInjurySetter(injuryn2,injury2)
            injury3 = milInjurySetter(injuryn3,injury3)
            
                hospitalNeed = injury.hospitalNeed
                if triage == 5
                    hospitalNeed = "H"
                    injury.hospitalNeed = "H"
                    injury2.hospitalNeed = "H"
                    injury3.hospitalNeed = "H"
                end

            if paramset.TriageLevel == "None"
                triage = 0
            end
            if consts.inputDBname=="vw24.sqlite"
                facility = convert(String,Patientslist[i,23])
                facilities = ["PEACH","APPLE","PEAR","WALNUT","PECAN","RAISIN","BANANA","BLACKBERRY","ORANGE","PLUM","MANGO","PINEAPPLE","CHESTNUT"]

                if facility in facilities
                
                    facility = "PoI"
                end
                if Patientslist[i,24] == 0
                    contam = false
                else
                    contam = true
                
                end
            end
            if consts.inputDBname=="ukr.sqlite"
                facility = "PoI"
                contam=false
            end

            patient = Patient(env, victag, age, is_treated, canwalk , lethally_injured,triage,x0,y0,iss,timeOfDeath,scheduledTime,SimParams,maxSimScore,bleed,(0,0),injury,injury2,injury3,hospitalNeed,consts,SS_array,distFromBlast,treatmentParams,0,facility,ip,contam)
           
           
            
            patient.hospitalneed=updateHospitalNeed(patient)
            if consts.inputDBname == "vw24.sqlite"
                if scheduledTime > 0.0
                    if victag == 23 
                        @yield timeout(env,scheduledTime)
                    end
                    if victag == 43 
                        @yield timeout(env,30.0)
                    end
                    if victag == 48 
                        @yield timeout(env,30.0)
                    end
                    if victag == 56 
                        @yield timeout(env,30.0)
                    end
                    if victag == 62
                        @yield timeout(env,30.0)
                    end
                    if victag == 69
                        @yield timeout(env,30.0)
                    end
                    if victag == 79
                        @yield timeout(env,30.0)
                    end
                    if victag == 85
                        @yield timeout(env,30.0)
                    end
                    if victag == 90
                        @yield timeout(env,30.0)
                    end
                    if victag == 101
                        @yield timeout(env,30.0)
                    end
                    if victag == 106
                        @yield timeout(env,30.0)
                    end
                    if victag == 112 
                        @yield timeout(env,30.0)
                    end
                    if victag == 134 
                        @yield timeout(env,30.0)
                    end
                    if victag == 141 
                        @yield timeout(env,30.0)
                    end
                    if victag == 157 
                        @yield timeout(env,30.0)
                    end
                    

                end
            end
            if consts.inputDBname == "ukr.sqlite"
                if scheduledTime > 0.0
                    if victag == 31
                        @yield timeout(env,scheduledTime)
                    end
                end
            end
            sscore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
                    
            if sscore <= 0.01 && patient.isdead == 0

                patient.triage = 9
                patient.isdead = 1
                VictimLogOutputDB(outputDB, patient.vicid,  now(env), 188,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,0,patient.SimParams[3],patient.bloodvolume)
                
            end
            if patient.triage !=9
                if consts.logg
                   
                    VictimLogOutputDB(outputDB, patient.vicid, now(env), 42,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,0,patient.SimParams[3],patient.bloodvolume)
                end
            end
       

            @yield request(patient.transport) #resource to let the victim function continue after transport

            # 
            if consts.inputDBname == "vw24.sqlite"
                ccp_data = Dict(
                    47.23234 => 1,  # CCP 1
                    47.22302 => 2,  # CCP 2
                    47.22367 => 3,  # CCP 3
                    47.23994 => 4,  # CCP 4
                    47.23822 => 5,  # CCP 5
                    47.23472 => 6,  # CCP 6
                    47.22055 => 7,  # CCP 7
                    47.2154  => 8,  # CCP 8
                    47.22353 => 9,  # CCP 9
                    47.220273 => 10 # CCP 10
                )

                default_idx = 1

                idx = get(ccp_data, Patientslist[i, 19], default_idx)
                

                @process victim(env, patient, outputDB, firefighters, unit,amst,mediclist, ambulist, medevac, transportqueue, transportqueueT3, hospitalQueue, paramset, treatment_data, 
                fmpUTM[idx], homeUTM[idx], ccpUTM[idx], fmp_LLA[idx], ccp_LLA[idx], timeCCPHome[idx], timeccpFMP[idx], timeFMPHome[idx], consts, inputDB, m, offroad)
            end
            if consts.inputDBname == "ukr.sqlite"
                if patient.vicid <= 8
                    @process victim(env, patient, outputDB, firefighters, unit,amst,mediclist, ambulist, medevac, transportqueue, transportqueueT3, hospitalQueue, paramset, treatment_data, 
                    fmpUTM[1], homeUTM[1], ccpUTM[1], fmp_LLA[1], ccp_LLA[1], timeCCPHome[1], timeccpFMP[1], timeFMPHome[1], consts, inputDB, m, offroad)
                else
                    @process victim(env, patient, outputDB, firefighters, unit,amst,mediclist, ambulist, medevac, transportqueue, transportqueueT3, hospitalQueue, paramset, treatment_data, 
                    fmpUTM[2], homeUTM[2], ccpUTM[2], fmp_LLA[2], ccp_LLA[2], timeCCPHome[2], timeccpFMP[2], timeFMPHome[2], consts, inputDB, m, offroad)
                end
            else
                @process victim(env, patient, outputDB, firefighters, unit,amst,mediclist, ambulist, medevac, transportqueue, transportqueueT3, hospitalQueue, paramset, treatment_data, 
                    fmpUTM[1], homeUTM[1], ccpUTM[1], fmp_LLA[1], ccp_LLA[1], timeCCPHome[1], timeccpFMP[1], timeFMPHome[1], consts, inputDB, m, offroad)
            end
        end
    end

    
    @process lastVictim(env, firefighters, mediclist, ambulist,medevac, paramset.Policy,paramset,outputDB,consts)
  
   
end
function infoHospitals(inputDB::SQLite.DB)

    hospitalInfo = unique(DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Hospital" )))
    return hospitalInfo
end
function computeMortality(outputDB::SQLite.DB,timeA::Float64,paramset::Parameterset)

    mortData = unique(DataFrame(DBInterface.execute( outputDB, "SELECT VictimID FROM VictimLog WHERE SimulationRun = $(paramset.simurun) AND Action = 188 AND Time <= $timeA" )))
    mort = size(mortData,1)

    return mort

end
function computeEvacuated(outputDB::SQLite.DB,timeA::Float64,paramset::Parameterset)

    mortData = unique(DataFrame(DBInterface.execute( outputDB, "SELECT VictimID FROM VictimLog WHERE SimulationRun = $(paramset.simurun) AND Action = 1002 AND Time <= $timeA" )))
    mort = size(mortData,1)
 
    return mort

end
function computeAction(outputDB::SQLite.DB,timeA::Float64,paramset::Parameterset,action::Int64)
    mortData = unique(DataFrame(DBInterface.execute( outputDB, "SELECT VictimID FROM VictimLog WHERE SimulationRun = $(paramset.simurun) AND Action = $action AND Time <= $timeA" )))
    mort = size(mortData,1)
 
    return mort
end
function infoVictims(outputDB::SQLite.DB,paramset::Parameterset,timeStart::Float64,timeEnd::Float64)
    victimPositions = (DataFrame(DBInterface.execute( outputDB, "SELECT * FROM VictimLog WHERE SimulationRun = $(paramset.simurun) AND Time <= $timeEnd AND Time >= $timeStart" )))
        
    VictimID = unique(victimPositions.VictimID)
    lastIndex =[findall(victimPositions.VictimID .== id )[end] for id in VictimID]
    victimPositions = victimPositions[lastIndex,:]
  
    return victimPositions
end
function infoTransport(outputDB::SQLite.DB,paramset::Parameterset,timeStart::Float64,timeEnd::Float64)
    victimPositions = (DataFrame(DBInterface.execute( outputDB, "SELECT * FROM VictimLog WHERE SimulationRun = $(paramset.simurun) AND (Action= 20 OR Action=17) AND Time <= $timeEnd AND Time >= $timeStart" )))
        
    VictimID = unique(victimPositions.VictimID)
    lastIndex =[findall(victimPositions.VictimID .== id )[end] for id in VictimID]
    victimPositions = victimPositions[lastIndex,:]
   
    return victimPositions
end

function addLayer(fm, victimPositions::DataFrame, transportPositions::DataFrame, quadrant::Int64, layerName::String)
    flm = pyimport("folium")
  
    fg = flm.FeatureGroup(name=layerName, show=false, overlay=false)
    lla_location = LLAfromUTM(quadrant, true, wgs84)
    
    # Add victim markers
    for j in 1:size(victimPositions, 1)
        patientUTM = Geodesy.UTM(victimPositions[j, :].patientX, victimPositions[j, :].patientY, 0.0)
        vicLLA = lla_location(patientUTM)
        
        if victimPositions[j, :].triage == 5 && victimPositions[j, :].Action != 23
            col = "purple"
            info = "ID: $(victimPositions[j, :].VictimID)\n<br>" *
                   "SS: $(victimPositions[j, :].SimScore)\n<br>" *
                   "triage: $(victimPositions[j, :].triage)\n<br>" *
                   "Action : $(victimPositions[j, :].Action)"
            flm.CircleMarker((vicLLA.lat - rand() * 0.0001, vicLLA.lon - rand() * 0.0001),
                popup=info,
                tooltip=info,
                radius=0.5,
                color=col,
                weight=10,
            ).add_to(fg)
        elseif victimPositions[j, :].triage == 3
            col = "green"
            info = "ID: $(victimPositions[j, :].VictimID)\n<br>" *
                   "SS: $(victimPositions[j, :].SimScore)\n<br>" *
                   "triage: $(victimPositions[j, :].triage)\n<br>" *
                   "Action : $(victimPositions[j, :].Action)"
            flm.CircleMarker((vicLLA.lat - rand() * 1e-4, vicLLA.lon - rand() * 1e-4),
                popup=info,
                tooltip=info,
                radius=0.5,
                color=col,
                weight=10,
            ).add_to(fg)
        elseif victimPositions[j, :].triage == 1
            col = "red"
            info = "ID: $(victimPositions[j, :].VictimID)\n<br>" *
                   "SS: $(victimPositions[j, :].SimScore)\n<br>" *
                   "triage: $(victimPositions[j, :].triage)\n<br>" *
                   "Action : $(victimPositions[j, :].Action)"
            flm.CircleMarker((vicLLA.lat - rand() * 1e-4, vicLLA.lon - rand() * 1e-4),
                popup=info,
                tooltip=info,
                radius=0.5,
                color=col,
                weight=10,
            ).add_to(fg)
        elseif victimPositions[j, :].triage == 2
            col = "orange"
            info = "ID: $(victimPositions[j, :].VictimID)\n<br>" *
                   "SS: $(victimPositions[j, :].SimScore)\n<br>" *
                   "triage: $(victimPositions[j, :].triage)\n<br>" *
                   "Action : $(victimPositions[j, :].Action)"
            flm.CircleMarker((vicLLA.lat - rand() * 1e-4, vicLLA.lon - rand() * 1e-4),
                popup=info,
                tooltip=info,
                radius=0.5,
                color=col,
                weight=10,
            ).add_to(fg)
        elseif victimPositions[j, :].triage == 9
            col = "black"
            info = "ID: $(victimPositions[j, :].VictimID)\n<br>" *
                   "SS: $(victimPositions[j, :].SimScore)\n<br>" *
                   "triage: $(victimPositions[j, :].triage)"
            flm.CircleMarker((vicLLA.lat, vicLLA.lon),
                popup=info,
                tooltip=info,
                radius=0.5,
                color=col,
                weight=10,
            ).add_to(fg)
        end
    end

    # Add transport markers
    for j in 1:size(transportPositions, 1)
        patientUTM = Geodesy.UTM(transportPositions[j, :].patientX, transportPositions[j, :].patientY, 0.0)
        vicLLA = lla_location(patientUTM)
        if transportPositions[j, :].Action == 17
            icon2 = flm.features.CustomIcon("ambuB.png", icon_size=(20, 20))
            info = "Ambulance dropping patient $(transportPositions[j, :].VictimID) at hospital $(transportPositions[j, :].hospid)\n<br>" *
               "lat: $(vicLLA.lat)\n<br>" *
               "lon: $(vicLLA.lon)"
            flm.Marker(
                (vicLLA.lat - rand() * 1e-4, vicLLA.lon - rand() * 1e-4),
                icon=icon2,
                tooltip=info
            ).add_to(fg)
        end
        if transportPositions[j, :].Action == 20
            icon2 = flm.features.CustomIcon("ambuB.png", icon_size=(20, 20))
            info = "Ambulance picking up patient $(transportPositions[j, :].VictimID) \n<br>" *
               "lat: $(vicLLA.lat)\n<br>" *
               "lon: $(vicLLA.lon)"
            flm.Marker(
                (vicLLA.lat - rand() * 1e-4, vicLLA.lon - rand() * 1e-4),
                icon=icon2,
                tooltip=info
            ).add_to(fg)
        end
       
    end

    if size(victimPositions, 1) > 0 || size(transportPositions, 1) > 0
        fg.add_to(fm)
    end
    
    flm.Marker((47.24625,18.12624),
        icon=flm.features.DivIcon(
        icon_size=(200,36),
        icon_anchor=(7,20),
        html="<div style='font-size: 18pt; color: red'>Situation at $layerName mins</div>"
        )  
    ).add_to(fg)
    
            
    return fm
end
function patientHosp(outputDB::SQLite.DB,id::Int64,paramset::Parameterset)
    
    patientResultsT1 = (DataFrame(DBInterface.execute( outputDB, "SELECT * FROM VictimLog WHERE SimulationRun = $(paramset.simurun) AND hospid = $id AND triage = 1 AND Action=17")))
    VictimID = unique(patientResultsT1.VictimID)
    lastIndex =[findall(patientResultsT1.VictimID .== id )[end] for id in VictimID]
    patientResultsT1 = patientResultsT1[lastIndex,:]
    nT1 = size(patientResultsT1,1)
    patientResultsT2 = (DataFrame(DBInterface.execute( outputDB, "SELECT * FROM VictimLog WHERE SimulationRun = $(paramset.simurun) AND hospid = $id AND triage = 2 AND Action=17")))
    VictimID = unique(patientResultsT2.VictimID)
    lastIndex =[findall(patientResultsT2.VictimID .== id )[end] for id in VictimID]
    patientResultsT2 = patientResultsT2[lastIndex,:]
    nT2 = size(patientResultsT2,1)
    return nT1,nT2
end

"""
    reportLocations(inputDB::SQLite.DB,outputDB::SQLite.DB,m::OpenStreetMapX.MapData,paramset::Parameterset,Action::Int64,routes::Vector{Vector{Int}},city::String,coords::Array,consts::Constants,timeStart::Float64 =-1.0,timeEnd::Float64=0.0,plotRoutes::Bool=false,plotHospitals::Bool=false)

    plots an html file, with locations of victims for different time intervals

"""
function reportLocations(inputDB::SQLite.DB,outputDB::SQLite.DB,m::OpenStreetMapX.MapData,paramset::Parameterset,Action::Int64,routes::Vector{Vector{Int}},city::String,coords::Array,consts::Constants,timeStart::Float64 =-1.0,timeEnd::Float64=0.0,plotRoutes::Bool=false,plotHospitals::Bool=false)
    
    flm = pyimport("folium");
    plugins = pyimport("folium.plugins")
   
    if timeStart >= 0.0

        overviewResults = (DataFrame(DBInterface.execute( outputDB, "SELECT * FROM Overview WHERE SimulationRun = $(paramset.simurun)")))
               
       #6 CSV.write("positions1.txt",victimPositions)
    end 
    fm = flm.Map()
        
    if timeStart >= 0.0 
        # feature_gp1 = flm.FeatureGroup(name="T1")
        # feature_gp2 = flm.FeatureGroup(name="T2")
        # feature_gp3 = flm.FeatureGroup(name="T3")
        # feature_gp4 = flm.FeatureGroup(name="dead")
        feature_gp5 = flm.FeatureGroup(name="BDE")
        feature_gp6 = flm.FeatureGroup(name="routes",show=false)
        feature_gp7 = flm.FeatureGroup(name="hospitals")
    
        feature_gp8 = flm.FeatureGroup(name="blast sites")
        
        #feature_gp9 = flm.FeatureGroup(name="Initial Situation")
        
        if plotRoutes #plot saved routes to file
            colrs = ["red","blue","green","purple","yellow","pink","orange","grey","black","white","brown"]
            for k in 1:length(routes)
                locs =[OpenStreetMapX.LLA(m.nodes[n],m.bounds) for n in routes[k]]
            
                info = "The route of ambulance $k\n<br>" *
                "Length: $(length(routes[k])) nodes\n<br>" *
                "From: $(routes[k][1]) $(round.((locs[1].lat, locs[1].lon),digits=4))\n<br>" *
                "To: $(routes[k][end]) $(round.((locs[end].lat, locs[end].lon),digits=4))"
                flm.PolyLine(
                    [(loc.lat, loc.lon) for loc in locs],
                    popup=info,
                    tooltip=info,
                    color=colrs[rand(1:length(colrs))]
            
                ).add_to(feature_gp6)
            end
           
        end
        if city != "kinshasa" && city !="waterloo"
            if city == "zapo"
                if plotHospitals
                    infoH = infoHospitals(inputDB)
                    count = 0
                    
                    for i in 1:size(infoH,1)
                        icon1 = flm.features.CustomIcon("role3.png", icon_size=(20,20))
                        info = "Name: $(infoH[i,:].Name)\n<br>" *
                        "id: $(infoH[i,:].ID)\n<br>" *
                        "capT1_A: $(infoH[i,:].CapT1_A_Medium)\n<br>" *
                        "T1,T2: $(patientHosp(outputDB,infoH[i,:].ID,paramset))"
                        t1,dum = patientHosp(outputDB,infoH[i,:].ID,paramset)
                        count = count + t1
                    
                        flm.Marker(
                            (infoH[i,:].LAT,infoH[i,:].LON),
                            icon=icon1,
                            tooltip=info
                        ).add_to(feature_gp7) 
                    end
                end  

            end
            if city == "cincu"
                if plotHospitals
                    infoH = infoHospitals(inputDB)
                    count = 0
                    
                    for i in 1:size(infoH,1)
                        if infoH[i,:].type == "R1"
                            icon1 = flm.features.CustomIcon("role1.png", icon_size=(20,20))
                            info = "Name: $(infoH[i,:].Name)\n<br>" *
                            "id: $(infoH[i,:].ID)\n<br>" *
                            "capT1_A: $(infoH[i,:].CapT1_A_Medium)\n<br>" *
                            "T1,T2: $(patientHosp(outputDB,infoH[i,:].ID,paramset))"
                            t1,dum = patientHosp(outputDB,infoH[i,:].ID,paramset)
                            count = count + t1
                        
                            flm.Marker(
                                (infoH[i,:].LAT,infoH[i,:].LON),
                                icon=icon1,
                                tooltip=info
                            ).add_to(feature_gp7) 

                          
                        end
                        if infoH[i,:].type == "R2B"
                            icon1 = flm.features.CustomIcon("role2.png", icon_size=(20,20))
                            info = "Name: $(infoH[i,:].Name)\n<br>" *
                            "id: $(infoH[i,:].ID)\n<br>" *
                            "capT1_A: $(infoH[i,:].CapT1_A_Medium)\n<br>" *
                            "T1,T2: $(patientHosp(outputDB,infoH[i,:].ID,paramset))"
                            t1,dum = patientHosp(outputDB,infoH[i,:].ID,paramset)
                            count = count + t1
                        
                            flm.Marker(
                                (infoH[i,:].LAT,infoH[i,:].LON),
                                icon=icon1,
                                tooltip=info
                            ).add_to(feature_gp7) 

                     
                        end
                        if infoH[i,:].type == "R3"
                            icon1 = flm.features.CustomIcon("role3.png", icon_size=(30,30))
                            info = "Name: $(infoH[i,:].Name)\n<br>" *
                            "id: $(infoH[i,:].ID)\n<br>" *
                            "capT1_A: $(infoH[i,:].CapT1_A_Medium)\n<br>" *
                            "T1,T2: $(patientHosp(outputDB,infoH[i,:].ID,paramset))"
                            t1,dum = patientHosp(outputDB,infoH[i,:].ID,paramset)
                            count = count + t1
                        
                            flm.Marker(
                                (infoH[i,:].LAT,infoH[i,:].LON),
                                icon=icon1,
                                tooltip=info
                            ).add_to(feature_gp7) 
                            
                            
                        end

                        for l in 1:length(coords[3])
                            icon2 = flm.features.CustomIcon("ccp1.png", icon_size=(20,20))
                        
                            
                            flm.Marker(
                                (coords[3][l].lat,coords[3][l].lon),
                                icon=icon2,
                            ).add_to(feature_gp7) 
                        end
                    end
                end  

            end
            if city == "vw24"

                
                flm.Polygon([(47.22355,18.17263),(47.22168,18.16803),(47.22399,18.16698),(47.22387,18.14962),(47.22692,18.14878),(47.23531,18.1606),(47.23071,18.16853)],  
                color="blue",
                weight=5,
                fill=true,
                fill_color="blue",
                fill_opacity=0.3,
                tooltip="BDE12"
                ).add_to(feature_gp5)
                flm.Polygon([(47.23517,18.16786),(47.2393,18.17025),(47.2408,18.17012),(47.24028,18.17628),(47.23074,18.18746),(47.22797,18.18235)],  
                color="blue",
                weight=5,
                fill=true,
                fill_color="blue",
                fill_opacity=0.3,
                tooltip="BDE11"
                ).add_to(feature_gp5)
                flm.Polygon([(47.22252,18.14901),(47.22926,18.13909),(47.23054,18.13353),(47.23116,18.12095),(47.22448,18.11956),(47.21204,18.13309),(47.21744,18.13968),(47.21953,18.13797),(47.22039,18.14444)],  
                color="blue",
                weight=5,
                fill=true,
                fill_color="blue",
                fill_opacity=0.3,
                tooltip="BDE13"
                ).add_to(feature_gp5)

                for k in 1:6
                    if k in 1:2
                    icon2 = flm.features.CustomIcon("ak.png", icon_size=(20,20))
                
                        info = "Event: TIC\n<br>" *
                        "lat: $(coords[1][k].lat)\n<br>" *
                        "lon: $(coords[1][k].lon)"
                
                        flm.Marker(
                            (coords[1][k].lat,coords[1][k].lon),
                            icon=icon2,
                            tooltip=info
                        ).add_to(feature_gp8) 
                    end
                    if k==3
                        icon2 = flm.features.CustomIcon("explo.png", icon_size=(30,30))
                    
                        info = "Event: IED\n<br>" *
                        "lat: $(coords[1][k].lat)\n<br>" *
                        "lon: $(coords[1][k].lon)"
                    
                        flm.Marker(
                            (coords[1][k].lat,coords[1][k].lon),
                            icon=icon2,
                            tooltip=info
                        ).add_to(feature_gp8) 
                    end
                    if k in 4:6
                        icon2 = flm.features.CustomIcon("c.png", icon_size=(20,20))
                   
                        info = "Event: Chem arty\n<br>" *
                        "lat: $(coords[1][k].lat)\n<br>" *
                        "lon: $(coords[1][k].lon)"
                    
                        flm.Marker(
                            (coords[1][k].lat,coords[1][k].lon),
                            icon=icon2,
                            tooltip=info
                        ).add_to(feature_gp8) 
                    end
                    icon2 = flm.features.CustomIcon("transport.png", icon_size=(20,20))
               
                
                    flm.Marker(
                        (47.23375,18.19443),
                        icon=icon2,
                    ).add_to(feature_gp7) 
                    for k in 1:10
                        icon2 = flm.features.CustomIcon("ccp1.png", icon_size=(20,20))
                       
                        
                        flm.Marker(
                            (coords[3][k].lat,coords[3][k].lon),
                            icon=icon2,
                        ).add_to(feature_gp7) 
                    end
                
                end
                # else #ruben demo
                #     for k in 1:3
                #         icon2 = flm.features.CustomIcon("C.png", icon_size=(20,20))
                
                #         info = "Event: GB attack\n<br>" *
                #         "lat: $(coords[1][k].lat)\n<br>" *
                #         "lon: $(coords[1][k].lon)"
                    
                #         flm.Marker(
                #             (coords[1][k].lat,coords[1][k].lon),
                #             icon=icon2,
                #             tooltip=info
                #         ).add_to(feature_gp8) 
                #     end

            end
       
        end
        if city == "waterloo"
                VictimsData = DataFrame(CSV.File("victimsArtillery.txt"; normalizenames=true))
                icon2 = flm.features.CustomIcon("crosshair.png", icon_size=(20,20)) #epicenter
               
                info = "Event: Arty strike\n<br>" *
                "lat: $(coords[1][1].lat)\n<br>" *
                "lon: $(coords[1][1].lon)"
            
                flm.Marker(
                    (coords[1][1].lat,coords[1][1].lon),
                    icon=icon2,
                    tooltip=info
                ).add_to(feature_gp8) 

        
                icon2 = flm.features.CustomIcon("ccp1.png", icon_size=(20,20))
         
                
                flm.Marker(
                    (coords[3][1].lat,coords[3][1].lon),
                    icon=icon2,
                ).add_to(feature_gp7) 
            

        end
    if city == "kinshasa"
            for k in 1:3
                icon2 = flm.features.CustomIcon("civ.png", icon_size=(20,20))
             
                
                flm.Marker(
                    (coords[3][k].lat,coords[3][k].lon),
                    icon=icon2,
                ).add_to(feature_gp7) 
            end

    end
        if plotHospitals && city == "vw24"
            infoH = infoHospitals(inputDB)
            count = 0
            if city != "vw24" && city != "cincu"
                for i in 1:size(infoH,1)
                    icon1 = flm.features.CustomIcon("role3.png", icon_size=(20,20))
                    info = "Name: $(infoH[i,:].Name)\n<br>" *
                    "id: $(infoH[i,:].ID)\n<br>" *
                    "capT1_A: $(infoH[i,:].CapT1_A_Medium)\n<br>" *
                    "T1,T2: $(patientHosp(outputDB,infoH[i,:].ID,paramset))"
                    t1,dum = patientHosp(outputDB,infoH[i,:].ID,paramset)
                    count = count + t1
                
                    flm.Marker(
                        (infoH[i,:].LAT,infoH[i,:].LON),
                        icon=icon1,
                        tooltip=info
                    ).add_to(feature_gp7) 
                end
            else
                for i in 1:size(infoH,1)
                    icon1 = flm.features.CustomIcon("role3.png", icon_size=(20,15))
                    if infoH[i,:].Name == "CAN"
                        icon1 = flm.features.CustomIcon("h40/ca.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "FRA"
                        icon1 = flm.features.CustomIcon("h40/fr.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "DNK"
                        icon1 = flm.features.CustomIcon("h40/dk.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "ARM"
                        icon1 = flm.features.CustomIcon("h40/am.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "GBR"
                        icon1 = flm.features.CustomIcon("h40/gb.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "ROU_MOL" || infoH[i,:].Name == "ROU_DSU" || infoH[i,:].Name == "ROU"
                        icon1 = flm.features.CustomIcon("h40/ro.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "ESP" || infoH[i,:].Name == "ESP_LTU"
                        icon1 = flm.features.CustomIcon("h40/es.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "KSV"
                        icon1 = flm.features.CustomIcon("h40/xk.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "NLD"
                        icon1 = flm.features.CustomIcon("h40/nl.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "BIH"
                        icon1 = flm.features.CustomIcon("h40/ba.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "USA" || infoH[i,:].Name == "USA_SVN"
                        icon1 = flm.features.CustomIcon("h40/us.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "BEL"
                        icon1 = flm.features.CustomIcon("h40/be.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "MNE"
                        icon1 = flm.features.CustomIcon("h40/me.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "BMTF_AZE"
                        icon1 = flm.features.CustomIcon("h40/bmtf1.png", icon_size=(30,30))
                    end
                    if infoH[i,:].Name == "HUN_RED"
                        icon1 = flm.features.CustomIcon("h40/hunred.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "LTU"
                        icon1 = flm.features.CustomIcon("h40/lt.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "EST_LTU" || infoH[i,:].Name == "EST"
                        icon1 = flm.features.CustomIcon("h40/ee.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "CZE"
                        icon1 = flm.features.CustomIcon("h40/cz.png", icon_size=(20,15))
                    end
                    if infoH[i,:].Name == "HUN"
                        icon1 = flm.features.CustomIcon("h40/hu.png", icon_size=(20,15))
                    end
                    info = "Name: $(infoH[i,:].Name)\n<br>" *
                    "id: $(infoH[i,:].ID)\n<br>" *
                    "capT1_A: $(infoH[i,:].CapT1_A_Medium)\n<br>" *
                    "T1,T2: $(patientHosp(outputDB,infoH[i,:].ID,paramset))"
                    t1,dum = patientHosp(outputDB,infoH[i,:].ID,paramset)
                    count = count + t1
                
                    flm.Marker(
                        (infoH[i,:].LAT,infoH[i,:].LON),
                        icon=icon1,
                        tooltip=info
                    ).add_to(feature_gp7) 
                end

            end
            
        end
    
        fm=addLayer(fm,infoVictims(outputDB,paramset,0.0,0.0),infoTransport(outputDB,paramset,0.0,0.0),consts.quadrant,"0") 
    

    # fm=addLayer(fm,infoVictims(outputDB,paramset,timeStart,timeEnd),quadrant,"Situation from $timeStart to $timeEnd")
        mort=[]
        for k in 1:120 #96 is for  480minutes
            fm=addLayer(fm,infoVictims(outputDB,paramset,0.0,k*5.0),infoTransport(outputDB,paramset,(k-1)*5.0,k*5.0),consts.quadrant,"$(k*5)")
            mortality=computeMortality(outputDB,k*5.0,paramset)
            push!(mort,[k*5 mortality])
            
        end
        
    end
    # inject html into the map html
    Policy=overviewResults[1,:].Policy
    victimsNumber=overviewResults[1,:].nrvictimcheck
    mortalityRun=overviewResults[1,:].nrdeathtotal
    mortality=computeMortality(outputDB,timeEnd,paramset)
    evacuated=computeEvacuated(outputDB,timeEnd,paramset)
    timesArray= range(0.0, stop = 600.0, step = 5.0)
    
    
    if city !="kinshasa"
        fm.get_root().html.add_child(flm.Element("""
        <div style="position: fixed; 
            top: 70px; left: 50px; width: 15%; height: 90%;background-color: #D3D3D3; 
            border:2px solid grey;z-index: 900;">
            <h5>SIMEDIS Info Box</h5>
            <h5>Policy : $Policy</h5>
            <h5>Number of victims : $victimsNumber</h5>
            <h5>Mortality at t=$timeEnd: $mortality</h5>
            <h5>Mortality at end of Run: $mortalityRun</h5>
            <h5>Victims admitted in H: $count</h5>
            <h5>Number of ambulances: $(paramset.nrAmbuFW)</h5>
            <ul style= margin-left: 1em">
                <li style="list-style-type: '&#128308';padding-left: 1em; padding-bottom: 0.25em">T1</li>
                <li style="list-style-type: '&#128992';padding-left: 1em; padding-bottom: 0.25em">T2</li>
                <li style="list-style-type: '&#128994';padding-left: 1em; padding-bottom: 0.25em">T3</li>
                <li style="list-style-type: '&#9899';padding-left: 1em; padding-bottom: 0.25em">Dead</li>
                <li style="list-style-image: url('ccp3.png');padding-left: 0.25em; padding-bottom: 0.25em">CCP</li>
            </ul>
            <button onclick="myFunction()">Plot</button>
            <p id="demo"></p>
            <!-- Create a canvas element for the chart -->
            <div style="position: relative; width: 100%; height: 15%; overflow: hidden;">
            <canvas id="myChart" style="position: absolute; left: 0; top: 0%;"></canvas>
            </div>
            <div style="position: relative; width: 100%; height: 15%; overflow: hidden;">
            <canvas id="myChart1" style="position: absolute; left: 0; top: 0%;"></canvas>
            </div>
            <div style="position: relative; width: 100%; height: 15%; overflow: hidden;">
            <canvas id="myChart2" style="position: absolute; left: 0; top: 0%;"></canvas>
            </div>
            <div style="position: relative; width: 100%; height: 15%; overflow: hidden;">
            <canvas id="myChart3" style="position: absolute; left: 0; top: 0%;"></canvas>
            </div>
            
            


            <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
            <script>
                // Set up the chart configuration
                var ctx = document.getElementById('myChart').getContext('2d');
                var ctx1 = document.getElementById('myChart1').getContext('2d');
                var ctx2 = document.getElementById('myChart2').getContext('2d');
                var ctx3 = document.getElementById('myChart3').getContext('2d');
                var chartConfig = {
                    type: 'line',
                    data: {
                        labels: [], // Will be populated with time values
                        datasets: [{
                            label: '# of dead vs time (mins)',
                            data: [], // Will be populated with mortality values
                            backgroundColor: 'blue'
                        }]
                    },
                    
                    options: {
                        responsive: true, // Set to false to disable responsiveness
                        maintainAspectRatio: false, // Set to false to disable aspect ratio preservation
                        scales: {
                            x: {
                                title: {
                                    display: true,
                                    text: 'Time (mins)'
                                }
                            },
                            y: {
                                title: {
                                    display: true,
                                    text: 'Mortality'
                                }
                            }
                        }
                    }
                };
                var chartConfigH = {
                    type: 'line',
                    data: {
                        labels: [], // Will be populated with time values
                        datasets: [{
                            label: '# arrived H vs time (mins)',
                            dataH: [], // Will be populated with mortality values
                            backgroundColor: 'green'
                        }]
                    },
                    
                    options: {
                        responsive: true, // Set to false to disable responsiveness
                        maintainAspectRatio: false, // Set to false to disable aspect ratio preservation
                        scales: {
                            x: {
                                title: {
                                    display: true,
                                    text: 'Time (mins)'
                                }
                            },
                            y: {
                                title: {
                                    display: true,
                                    text: 'Arrived H'
                                }
                            }
                        }
                    }
                };
                var chartConfigFMP = {
                    type: 'line',
                    data: {
                        labels: [], // Will be populated with time values
                        datasets: [{
                            label: '# arrived FMP vs time (mins)',
                            dataFMP: [], // Will be populated with mortality values
                            backgroundColor: 'red'
                        }]
                    },
                    
                    options: {
                        responsive: true, // Set to false to disable responsiveness
                        maintainAspectRatio: false, // Set to false to disable aspect ratio preservation
                        scales: {
                            x: {
                                title: {
                                    display: true,
                                    text: 'Time (mins)'
                                }
                            },
                            y: {
                                title: {
                                    display: true,
                                    text: 'Arrived FMP'
                                }
                            }
                        }
                    }
                };
                var chartConfigCCP = {
                    type: 'line',
                    data: {
                        labels: [], // Will be populated with time values
                        datasets: [{
                            label: '# Arrived CCP vs time (mins)',
                            dataCCP: [], // Will be populated with victim values
                            backgroundColor: 'orange'
                        }]
                    },
                    
                    options: {
                        responsive: true, // Set to false to disable responsiveness
                        maintainAspectRatio: false, // Set to false to disable aspect ratio preservation
                        scales: {
                            x: {
                                title: {
                                    display: true,
                                    text: 'Time (mins)'
                                }
                            },
                            y: {
                                title: {
                                    display: true,
                                    text: 'Arrived CCP'
                                }
                            }
                        }
                    }
                };

                // Create the chart
                var myChart = new Chart(ctx, chartConfig);
                var myChart1 = new Chart(ctx1, chartConfigH);
                var myChart2 = new Chart(ctx2, chartConfigFMP);
                var myChart3 = new Chart(ctx3, chartConfigCCP);
               
                function myFunction() {
                    // Define your data directly
                    var data = [];
                    var times = [0.0, 5.0, 10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,105.0,110.0,115.0];; // Define your time values here
                    var mortalitys = [$(computeMortality(outputDB, 0.0, paramset)),$(computeMortality(outputDB, 5.0, paramset)),$(computeMortality(outputDB, 10.0, paramset)),$(computeMortality(outputDB, 15.0, paramset)),$(computeMortality(outputDB, 20.0, paramset)),$(computeMortality(outputDB, 25.0, paramset)),$(computeMortality(outputDB, 30.0, paramset)),$(computeMortality(outputDB, 35.0, paramset)),$(computeMortality(outputDB, 40.0, paramset)),$(computeMortality(outputDB, 45.0, paramset)),$(computeMortality(outputDB, 50.0, paramset)),$(computeMortality(outputDB, 55.0, paramset)),$(computeMortality(outputDB, 60.0, paramset)),$(computeMortality(outputDB, 65.0, paramset)),$(computeMortality(outputDB, 70.0, paramset)),$(computeMortality(outputDB, 75.0, paramset)),$(computeMortality(outputDB, 80.0, paramset)),$(computeMortality(outputDB, 85.0, paramset)),$(computeMortality(outputDB, 90.0, paramset)),$(computeMortality(outputDB, 95.0, paramset)),$(computeMortality(outputDB, 100.0, paramset)),$(computeMortality(outputDB, 105.0, paramset)),$(computeMortality(outputDB, 110.0, paramset)),$(computeMortality(outputDB, 115.0, paramset))]


                    for (var i = 0; i < times.length; i++) {
                        var time = times[i]
                        var mortality = mortalitys[i];
                        data.push({ time: time, mortality: mortality });
                    }
                    var dataH = [];
                   

                    var timesH = [0.0, 5.0, 10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,105.0,110.0,115.0]; // Define your time values here
                    var Harrivals = [$(computeAction(outputDB, 0.0, paramset,17)),$(computeAction(outputDB, 5.0, paramset,17)),$(computeAction(outputDB, 10.0, paramset,17)),$(computeAction(outputDB, 15.0, paramset,17)),$(computeAction(outputDB, 20.0, paramset,17)),$(computeAction(outputDB, 25.0, paramset,17)),$(computeAction(outputDB, 30.0, paramset,17)),$(computeAction(outputDB, 35.0, paramset,17)),$(computeAction(outputDB, 40.0, paramset,17)),$(computeAction(outputDB, 45.0, paramset,17)),$(computeAction(outputDB, 50.0, paramset,17)),$(computeAction(outputDB, 55.0, paramset,17)),$(computeAction(outputDB,60.0, paramset,17)),$(computeAction(outputDB, 65.0, paramset,17)),$(computeAction(outputDB, 70.0, paramset,17)),$(computeAction(outputDB, 75.0, paramset,17)),$(computeAction(outputDB, 80.0, paramset,17)),$(computeAction(outputDB, 85.0, paramset,17)),$(computeAction(outputDB, 90.0, paramset,17)),$(computeAction(outputDB, 95.0, paramset,17)),$(computeAction(outputDB, 100.0, paramset,17)),$(computeAction(outputDB, 105.0, paramset,17)),$(computeAction(outputDB, 110.0, paramset,17)),$(computeAction(outputDB, 115.0, paramset,17))]


                    for (var i = 0; i < timesH.length; i++) {
                        var timeH = timesH[i]
                        var Harrival = Harrivals[i];
                        dataH.push({ timeH: timeH, Harrival: Harrival });
                    }
                    var dataFMP = [];
                    var timesFMP = [0.0, 5.0, 10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,105.0,110.0,115.0]; // Define your time values here
                    var FMParrivals = [$(computeAction(outputDB, 0.0, paramset,38)),$(computeAction(outputDB, 5.0, paramset,38)),$(computeAction(outputDB, 10.0, paramset,38)),$(computeAction(outputDB, 15.0, paramset,38)),$(computeAction(outputDB, 20.0, paramset,38)),$(computeAction(outputDB, 25.0, paramset,38)),$(computeAction(outputDB, 30.0, paramset,38)),$(computeAction(outputDB, 35.0, paramset,38)),$(computeAction(outputDB, 40.0, paramset,38)),$(computeAction(outputDB, 45.0, paramset,38)),$(computeAction(outputDB, 50.0, paramset,38)),$(computeAction(outputDB, 55.0, paramset,38)),$(computeAction(outputDB,60.0, paramset,38)),$(computeAction(outputDB, 65.0, paramset,38)),$(computeAction(outputDB, 70.0, paramset,38)),$(computeAction(outputDB, 75.0, paramset,38)),$(computeAction(outputDB, 80.0, paramset,38)),$(computeAction(outputDB, 85.0, paramset,38)),$(computeAction(outputDB, 90.0, paramset,38)),$(computeAction(outputDB, 95.0, paramset,38)),$(computeAction(outputDB, 100.0, paramset,38)),$(computeAction(outputDB, 105.0, paramset,38)),$(computeAction(outputDB, 110.0, paramset,38)),$(computeAction(outputDB, 115.0, paramset,38))]


                    for (var i = 0; i < timesFMP.length; i++) {
                        var timeFMP = timesFMP[i]
                        var FMParrival = FMParrivals[i];
                        dataFMP.push({ timeFMP: timeFMP, FMParrival: FMParrival });
                    }
                    var dataCCP = [];
                    var timesCCP = [0.0, 5.0, 10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,105.0,110.0,115.0]; // Define your time values here
                    var CCParrivals = [$(computeAction(outputDB, 0.0, paramset,1)),$(computeAction(outputDB, 5.0, paramset,1)),$(computeAction(outputDB, 10.0, paramset,1)),$(computeAction(outputDB, 15.0, paramset,1)),$(computeAction(outputDB, 20.0, paramset,1)),$(computeAction(outputDB, 25.0, paramset,1)),$(computeAction(outputDB, 30.0, paramset,1)),$(computeAction(outputDB, 35.0, paramset,1)),$(computeAction(outputDB, 40.0, paramset,1)),$(computeAction(outputDB, 45.0, paramset,1)),$(computeAction(outputDB, 50.0, paramset,1)),$(computeAction(outputDB, 55.0, paramset,1)),$(computeAction(outputDB,60.0, paramset,1)),$(computeAction(outputDB, 65.0, paramset,1)),$(computeAction(outputDB, 70.0, paramset,1)),$(computeAction(outputDB, 75.0, paramset,1)),$(computeAction(outputDB, 80.0, paramset,1)),$(computeAction(outputDB, 85.0, paramset,1)),$(computeAction(outputDB, 90.0, paramset,1)),$(computeAction(outputDB, 95.0, paramset,1)),$(computeAction(outputDB, 100.0, paramset,1)),$(computeAction(outputDB, 105.0, paramset,1)),$(computeAction(outputDB, 110.0, paramset,1)),$(computeAction(outputDB, 115.0, paramset,1))]


                    for (var i = 0; i < timesCCP.length; i++) {
                        var timeCCP = timesCCP[i]
                        var CCParrival = CCParrivals[i];
                        dataCCP.push({ timeCCP: timeCCP, CCParrival: CCParrival });
                    }


                    // Update the chart with the defined data
                    updateChart(data,dataH,dataFMP,dataCCP);
                }

                function updateChart(data,dataH,dataFMP,dataCCP) {
                    var times = data.map(entry => entry.time);
                    var mortality = data.map(entry => entry.mortality);
                    var timesH = dataH.map(entry => entry.timeH);
                    var Harrival = dataH.map(entry => entry.Harrival);
                    var timesFMP = dataFMP.map(entry => entry.timeFMP);
                    var FMParrival = dataFMP.map(entry => entry.FMParrival);
                    var timesCCP = dataCCP.map(entry => entry.timeCCP);
                    var CCParrival = dataCCP.map(entry => entry.CCParrival);

                    // Update the chart data and labels
                    myChart.data.labels = times;
                    myChart.data.datasets[0].data = mortality;
                    myChart.update();
                    // Create a new dataset for myChart1 and update it with dataH
                    var datasetH = {
                        label: '# arrived H vs time (mins)',
                        data: Harrival,
                        backgroundColor: 'green'
                    };
                    myChart1.data.labels = timesH;
                    myChart1.data.datasets = [datasetH]; // Replace the datasets array
                    myChart1.update();

                    // Create a new dataset for myChart2 and update it with dataFMP
                    var datasetFMP = {
                        label: '# arrived FMP vs time (mins)',
                        data: FMParrival,
                        backgroundColor: 'red'
                    };
                    myChart2.data.labels = timesFMP;
                    myChart2.data.datasets = [datasetFMP]; // Replace the datasets array
                    myChart2.update();

                    var datasetCCP = {
                        label: '# arrived CCP vs time (mins)',
                        data: CCParrival,
                        backgroundColor: 'orange'
                    };
                    myChart3.data.labels = timesCCP;
                    myChart3.data.datasets = [datasetCCP]; // Replace the datasets array
                    myChart3.update();

                    
                    // Update the "demo" element
                    document.getElementById("demo").innerHTML = "Charts updated.";
                }
            </script>
        </div>
        <div style="position: fixed; 
        top: 85%; left: 3%; width: 7%; height: 7%;z-index: 900;">
             <img src="defence1.png" width=50%><img src="khid.png" width=50%>
            
        </div>
        
        
        """))
    end
    if city == "kinshasa"
        fm.get_root().html.add_child(flm.Element("""
        <div style="position: fixed; 
            top: 70px; left: 70px; width: 300px; height: 1000px;background-color: #D3D3D3; 
            border:2px solid grey;z-index: 900;">
            <h5>SIMEDIS Info Box</h5>
            <h5>Policy : $Policy</h5>
            <h5>Number of victims : $victimsNumber</h5>
            <h5>Evacuated at t=$timeEnd: $evacuated</h5>
            <h5>Number of boats : $(consts.nrBoats) / Boat capacity: $(consts.BoatCapacity)</h5>
            <ul style= margin-left: 1em">
                <li style="list-style-type: '&#128308';padding-left: 1em; padding-bottom: 0.25em">T1</li>
                <li style="list-style-type: '&#128992';padding-left: 1em; padding-bottom: 0.25em">T2</li>
                <li style="list-style-type: '&#128994';padding-left: 1em; padding-bottom: 0.25em">T3</li>
                <li style="list-style-type: '&#9899';padding-left: 1em; padding-bottom: 0.25em">Dead</li>
                <li style="list-style-type: '&#128995';padding-left: 1em; padding-bottom: 0.25em">refugee</li>
            </ul>
            <button onclick="myFunction()">Plot</button>
            <p id="demo"></p>
            <!-- Create a canvas element for the chart -->
            <div style="position: relative; width: 100%; height: 200px; overflow: hidden;">
            <canvas id="myChart" style="position: absolute; left: 0; top: 0;"></canvas>
            </div>
        

            <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
            <script>
                // Set up the chart configuration
                var ctx = document.getElementById('myChart').getContext('2d');
                var chartConfig = {
                    type: 'line',
                    data: {
                        labels: [], // Will be populated with time values
                        datasets: [{
                            label: '# of evacuated vs time (mins)',
                            data: [], // Will be populated with mortality values
                            backgroundColor: 'blue'
                        }]
                    },
                    options: {
                        responsive: true, // Set to false to disable responsiveness
                        maintainAspectRatio: false, // Set to false to disable aspect ratio preservation
                        scales: {
                            x: {
                                title: {
                                    display: true,
                                    text: 'Time (mins)'
                                }
                            },
                            y: {
                                title: {
                                    display: true,
                                    text: '# Evacuated'
                                }
                            }
                        }
                    }
                };

                // Create the chart
                var myChart = new Chart(ctx, chartConfig);

                function myFunction() {
                    // Define your data directly
                    var data = [];
                    var times = [0.0,3.0, 5.0, 10.0,12.0,15.0,20.0,25.0,30.0,35.0,40.0]; // Define your time values here
                    var mortalitys = [$(computeEvacuated(outputDB, 0.0, paramset)),$(computeEvacuated(outputDB, 3.0, paramset)),$(computeEvacuated(outputDB, 5.0, paramset)),$(computeEvacuated(outputDB, 10.0, paramset)),$(computeEvacuated(outputDB, 12.0, paramset)),$(computeEvacuated(outputDB, 15.0, paramset)),$(computeEvacuated(outputDB, 20.0, paramset)),$(computeEvacuated(outputDB, 25.0, paramset)),$(computeEvacuated(outputDB, 30.0, paramset)),$(computeEvacuated(outputDB, 35.0, paramset)),$(computeEvacuated(outputDB, 40.0, paramset)),]


                    for (var i = 0; i < times.length; i++) {
                        var time = times[i]
                        var mortality = mortalitys[i];
                        data.push({ time: time, mortality: mortality });
                    }

                    // Update the chart with the defined data
                    updateChart(data);
                }

                function updateChart(data) {
                    var times = data.map(entry => entry.time);
                    var mortality = data.map(entry => entry.mortality);

                    // Update the chart data and labels
                    myChart.data.labels = times;
                    myChart.data.datasets[0].data = mortality;
                    myChart.update();

                    // Update the "demo" element
                    document.getElementById("demo").innerHTML = "Chart updated.";
                }
            </script>
        </div>
        
        
        """))


    end
    # add tile layers
    flm.TileLayer("openstreetmap").add_to(fm)
    flm.TileLayer(tiles = "https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}",
    attr = "Google",
    name = "Google Satellite",
    overlay = true,
    control = true
    ).add_to(fm)

    
    flm.TileLayer(
        tiles = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        attr = "i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community",
        name = "Esri Satellite",
        overlay = true,
        control = true
    ).add_to(fm)
  
    plugins.MousePosition().add_to(fm)
    plugins.Geocoder(position="bottomleft").add_to(fm)

 
   # flm.TileLayer("stamenterrain", attr="stamenterrain").add_to(fm)

    
    # add layers control over the map
  
   
    feature_gp5.add_to(fm)
    feature_gp6.add_to(fm)
    feature_gp7.add_to(fm)
    feature_gp8.add_to(fm)
   
   # feature_gp9.add_to(fm)


    flm.LayerControl().add_to(fm)
    MAP_BOUNDS = [(m.bounds.min_y,m.bounds.min_x),(m.bounds.max_y,m.bounds.max_x)]
    flm.Rectangle(MAP_BOUNDS, color="black", weight=6).add_to(fm)
    fm
    fm.fit_bounds(MAP_BOUNDS)
    
    # text = "this is where you can read all the info you need"

    # iframe = flm.IFrame(text, width=700, height=450)
    # popup = flm.Popup(iframe, max_width=3000)
    
    # Text = flm.Marker((-4.311607617931119, 15.295258066453652), popup=popup,
    #                     icon=flm.Icon(icon_color="green"))
    # fm.add_child(Text)
   
    #   if Action == 17 || Action == 188 && timeStart == 0.0
        
    #     fm.save("Situation$(city)$(Action)_$(paramset.Policy)_$(paramset.HospitalDistribution).html")
    #   end
   
  
  if timeStart >= 0.0
    fm.save("Situation-$(city)$(paramset.simurun)$(paramset.HospitalDistribution)$(paramset.Policy)Tourniquet$(paramset.Tourniquet)Ambulance$(paramset.nrAmbuFW)_$plotRoutes.html")
  end
 
  

end

@resumable function spawnVictims(env ::Environment, event:: Event,SS_array::Array,treatment_data::Array,location::UTM,outputDB ::SQLite.DB, firefighters ::Resource, boatList::BoatList,mediclist ::MedicList, ambulist ::AmbuList, medevac::Medevac, transportqueue ::Array,transportqueue1::Array, transportqueueT3 ::Array, hospitalQueue ::Array, paramset ::Parameterset, consts::Constants,inputDB,m,fmp_LLA,ccp_LLA,homeUTM,ccpUTM,timeCCPHome,timeccpFMP,timeFMPHome,fmpUTM,threat,offroad)
   #this function spawns numberofVictims victims at $location()  
   unit = DeconUnit(env,20,1000,40,100,1000,300,1000,100)
   amst = Amst(env,100,40,false)
   victims_list0 = []
   victims_list1 = []
   victims_list2 = []
   victims_list=[]
   d = Normal(30, 10)
   ageArray = rand(d, 1000) #generate 1000 points from an age pyramid according to a normal law with 40 as the mean age and 20 as the standard deviation
   if consts.creator 
       if event.cbrn==5
        if threat ==0
            push!(victims_list0, ["victag" "age" "triage" "threat" "hitX" "hitY" "x0"  "y0" "distFromBlast" "ISS"  "bGompertz"  "cGompertz"   "gamma" "toD" "patient.mobile" "bleed" "viclat"  "viclon"  "hitlat"    "hitlon" "scheduledTime"])
        end
        if threat==1
            push!(victims_list1, ["victag" "age" "triage" "threat" "hitX" "hitY" "x0"  "y0" "distFromBlast" "ISS"  "bGompertz"  "cGompertz"   "gamma" "toD" "patient.mobile" "bleed" "viclat"  "viclon"  "hitlat"    "hitlon" "scheduledTime"])
        end
        if threat==2
            push!(victims_list2, ["victag" "age" "triage" "threat" "hitX" "hitY" "x0"  "y0" "distFromBlast" "ISS"  "bGompertz"  "cGompertz"   "gamma" "toD" "patient.mobile" "bleed" "viclat"  "viclon"  "hitlat"    "hitlon" "scheduledTime"])
        end
       end
       if event.cbrn==4
        
        push!(victims_list, ["victag" "age" "triage" "threat" "hitX" "hitY" "x0"  "y0" "distFromBlast" "ISS"  "bGompertz"  "cGompertz"   "gamma" "toD" "patient.mobile" "bleed" "brad_Gompertz" "crad_Gompertz" "gammarad_Gompertz" "tDeath_rad" "dose" "scheduledTime"])
       end
       if event.cbrn==0 #refugeeEvac
        if threat ==0
            push!(victims_list0, ["victag" "age" "triage" "threat" "hitX" "hitY" "x0"  "y0" "ISS"  "bGompertz"  "cGompertz"   "gamma" "toD" "patient.mobile" "bleed" "viclat"  "viclon"  "hitlat"    "hitlon" "scheduledTime"])
        end
        if threat==1
            push!(victims_list1, ["victag" "age" "triage" "threat" "hitX" "hitY" "x0"  "y0" "ISS"  "bGompertz"  "cGompertz"   "gamma" "toD" "patient.mobile" "bleed" "viclat"  "viclon"  "hitlat"    "hitlon" "scheduledTime"])
        end
        if threat==2
            push!(victims_list2, ["victag" "age" "triage" "threat" "hitX" "hitY" "x0"  "y0" "ISS"  "bGompertz"  "cGompertz"   "gamma" "toD" "patient.mobile" "bleed" "viclat"  "viclon"  "hitlat"    "hitlon" "scheduledTime"])
        end
       end

   end
   if (event.cbrn==5) || (event.cbrn==0)
    victimPos = generateVictimPosition(env, event,location,consts)
   
   end
   if (event.cbrn==4)
    victimPos = generateVictimPositionRad(env, event,location,consts)
   end
   
  

   for i in 1:event.casualtyEstimate #generate numberofVictims victims
       
       victag = i+event.casualtyEstimate*threat
       is_treated = 0
       lethally_injured = 0
       agef = abs(ageArray[i])
       if event.cbrn==0
        age=max(0,floor(Int64,agef))
       else
        age=max(18,floor(Int64,agef))
       end
       if event.cbrn==0
        age=max(0,floor(Int64,agef))
       else
        age=max(18,floor(Int64,agef))
       end
       x0 = victimPos[i,1]
       y0 = victimPos[i,2]
       xLat = victimPos[i,5]
       xLon = victimPos[i,6]
       hitLat = victimPos[i,7]
       hitLon = victimPos[i,8]
             
       xLat = victimPos[i,5]
       xLon = victimPos[i,6]
       hitLat = victimPos[i,7]
       hitLon = victimPos[i,8]
             
       armor = 100.0
       maxSimScore = 20.0
       injury =Injury(1,"","","",1,[""],"",6.0,"",0.0)
       injury2 =Injury(1,"","","",1,[""],"",6.0,"",0.0)
       injury3 =Injury(1,"","","",1,[""],"",6.0,"",0.0)

       injuryn = 0
       injuryn2= 0
       injuryn3 = 0
       
       injury = milInjurySetter(injuryn,injury)
       injury2 = milInjurySetter(injuryn2,injury2)
       injury3 = milInjurySetter(injuryn3,injury3)
       global issSum = 0
       canwalk = 1
       bleed = 0
       #For each blast
        if(event.cbrn==5) 
            for j in 1:event.casualtyEstimate
                
                if j % event.casualtyRate == 0
                 
                 
                    hitloc = (victimPos[j,3],victimPos[j,4])
                    
                    distFromBlast = sqrt((hitloc[1]-x0)^2+(hitloc[2]-y0)^2)
                                
                #Initial victim Location close from the threat
                    lethal=rand()

                
                    canwalk1,probIncap,incapRand=walkable(distFromBlast,event.aoe)

                    if canwalk1 == 0 #if once unable to walk always unable to
                 
                 
                        canwalk = 0
                                                            
                    end
                    
                    penetrationThreshold = shrapnelPenetration(distFromBlast,10000)
                    
                    penRand = rand()
                    penetrationThreshold = penetrationThreshold/2.0 #supposed mitigation over all body
                    if penRand < penetrationThreshold && bleed == 0
                        
                        hi = "hit"
                        hip = 1
                        bleed = 1
                        lethally_injured=1 #hemorrage is deadly
                        
                    else
                        
                        hi="dodge"
                        hip = 0
                        
                
                    end
                    if canwalk == 1
                        incap = 0
                    end
                    if canwalk == 0 && bleed == 1
                        incap = 1
                    end
                    if distFromBlast < 15 || (hip == 1 && incap ==1) 
                        lethally_injured=1
                        
                    end
                    
                    if distFromBlast >= 15 && distFromBlast < 24
                        if lethal <= 0.3 || (hip == 1 && incap == 1) 
                            lethally_injured=1
                        else
                            
                        end
                                
                    end
                    if distFromBlast >= 24 && distFromBlast< 40

                        if lethal <= 0.05 || (hip == 1 && incap ==1) 
                            lethally_injured=1
                        else
                            
                        end
                    end
                    if distFromBlast >= 40 && distFromBlast < 100

                        if lethal <= 0.01 || (hip == 1 && incap ==1) 
                            lethally_injured=1
                        
                        else
                            
                        end
                    end
                    if distFromBlast >= 100 

                        if lethal <= 0.001 || (hip == 1 && incap ==1) 
                            lethally_injured=1
                        else
                            
                        end
                    end
                    if event.name == "artillery155mm"
                        iss = trunc(Int64,min(75.0,3*75.0/(distFromBlast))) #artillery
                    end
                    if event.name == "droneStrike"
                        iss = trunc(Int64,min(75.0,2*75.0/(distFromBlast))) #drone
                    end
                    if event.name == "FPVDrone"
                        iss = trunc(Int64,min(75.0,1*75.0/(distFromBlast))) #drone
                    end

                    if event.name == "Chem"
                        iss = trunc(Int64,min(75.0,(0.66)*75.0/(distFromBlast))) #drone
                    end
                    
                
                    global issSum = issSum + iss
                
                end
                
            end
            iss = trunc(Int64,min(75.0,issSum))
    
         
           
            if iss > 0
                timeOfDeath = 43500*iss^(-1.95)
            else 
                timeOfDeath = 1.0e6
                iss = 0
            end
            ip = 7 
            SimParams=[]
            SimParams=[]
            triage = getTriage(iss,ip)
            treatmentParams=[0.0,0.0,0.0,0.0,0.0]
            treatmentParams=[0.0,0.0,0.0,0.0,0.0]
            
            hospitalNeed = injury.hospitalNeed
            armor = 100.0
            #scheduledTime = 0.0
            facility="PoI"
            patient = Patient(env, victag, age, is_treated, canwalk , lethally_injured,triage,x0,y0,iss,timeOfDeath,event.scheduledTime,SimParams,maxSimScore,bleed,(0,0),injury,injury2,injury3,hospitalNeed,consts,SS_array,distFromBlast,treatmentParams,threat,facility,ip)
            patient.hospitalneed=updateHospitalNeed(patient)
         

           
            if (lethally_injured == 1)
                patient.SimParams[3]=1.0
            end
            
            realToD = 20000.0

            for i in 1:2000
                test = checkHealth(env,patient,convert(Float64,i),0.01)
                if test == 1
                    realToD = convert(Float64,i)
                    
                    break
                end
            end
            
            if consts.creator
               # push!(victims,patient)
               if threat == 0
                push!(victims_list0, [victag floor(Int64,age) triage event.name round(hitloc[1],digits=2) round(hitloc[2],digits=2)  round(patient.x0,digits=2)  round(patient.y0,digits=2) round(distFromBlast,digits=2) floor(Int64,iss)  round(patient.SimParams[1],digits=2)  round(patient.SimParams[2],digits=2)    patient.SimParams[3] round(patient.SimParams[4],digits=2) patient.mobile bleed xLat xLon hitLat hitLon event.scheduledTime])
               end
               if threat == 1
                push!(victims_list1, [victag floor(Int64,age) triage event.name round(hitloc[1],digits=2) round(hitloc[2],digits=2)  round(patient.x0,digits=2)  round(patient.y0,digits=2) round(distFromBlast,digits=2) floor(Int64,iss)  round(patient.SimParams[1],digits=2)  round(patient.SimParams[2],digits=2)    patient.SimParams[3] round(patient.SimParams[4],digits=2) patient.mobile bleed xLat xLon hitLat hitLon event.scheduledTime])
               end
               if threat == 2
                push!(victims_list2, [victag floor(Int64,age) triage event.name round(hitloc[1],digits=2) round(hitloc[2],digits=2)  round(patient.x0,digits=2)  round(patient.y0,digits=2) round(distFromBlast,digits=2) floor(Int64,iss)  round(patient.SimParams[1],digits=2)  round(patient.SimParams[2],digits=2)    patient.SimParams[3] round(patient.SimParams[4],digits=2) patient.mobile bleed xLat xLon hitLat hitLon event.scheduledTime])
               end
            end
            if event.scheduledTime >= 0.0 && victag == 31
                @yield timeout(env,event.scheduledTime)
            end

            sscore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
          
            @yield request(patient.transport) #resource to let the victim function continue after transport
            @process victim(env, patient, outputDB, firefighters,unit,amst, mediclist, ambulist,medevac, transportqueue, transportqueueT3, hospitalQueue, paramset,treatment_data,fmpUTM,homeUTM,ccpUTM,fmp_LLA,ccp_LLA,timeCCPHome,timeccpFMP,timeFMPHome,consts,inputDB,m,offroad)
        end 
        if (event.cbrn==0)
            for j in 1:event.casualtyEstimate 
                
                    hitloc =(0.0,0.0)
                    
                    distFromBlast = 0.0 #irrelevant
                    canwalk = 1 #byDefault     
                    bleed=0
                    
                    lethally_injured = 0
                    iss = 0
            end
 
            if iss > 0
                timeOfDeath = 43500*iss^(-1.95)
            else 
                timeOfDeath = 1.0e6
                iss = 0
            end
            ip = 7 
            SimParams=[]
            triage = getTriage(iss,ip)
            treatmentParams=[0.0,0.0,0.0,0.0,0.0]
            
            hospitalNeed = injury.hospitalNeed
            armor = 100.0
            scheduledTime = 0.0
            facility="PoI"
            patient = Patient(env, victag, age, is_treated, canwalk , lethally_injured,triage,x0,y0,iss,timeOfDeath,scheduledTime,SimParams,maxSimScore,bleed,(0,0),injury,injury2,injury3,hospitalNeed,consts,SS_array,distFromBlast,treatmentParams,threat,facility,ip)
            patient.hospitalneed=updateHospitalNeed(patient)
            if (lethally_injured == 1)
                patient.SimParams[3]=1.0
            end
            
            realToD = 20000.0

            for i in 1:2000
                test = checkHealth(env,patient,convert(Float64,i),0.01)
                if test == 1
                    realToD = convert(Float64,i)
                    
                    break
                end
            end
            
            if consts.creator
               # push!(victims,patient)
               if threat == 0
                push!(victims_list0, [victag floor(Int64,age) triage event.name round(hitloc[1],digits=2) round(hitloc[2],digits=2)  round(patient.x0,digits=2)  round(patient.y0,digits=2) round(distFromBlast,digits=2) floor(Int64,iss)  round(patient.SimParams[1],digits=2)  round(patient.SimParams[2],digits=2)    patient.SimParams[3] round(patient.SimParams[4],digits=2) patient.mobile bleed xLat xLon hitLat hitLon])
               end
               if threat == 1
                push!(victims_list1, [victag floor(Int64,age) triage event.name round(hitloc[1],digits=2) round(hitloc[2],digits=2)  round(patient.x0,digits=2)  round(patient.y0,digits=2) floor(Int64,iss) round(distFromBlast,digits=2)   round(patient.SimParams[1],digits=2)  round(patient.SimParams[2],digits=2)    patient.SimParams[3] round(patient.SimParams[4],digits=2) patient.mobile bleed xLat xLon hitLat hitLon])
               end
               if threat == 2
                push!(victims_list2, [victag floor(Int64,age) triage event.name round(hitloc[1],digits=2) round(hitloc[2],digits=2)  round(patient.x0,digits=2)  round(patient.y0,digits=2) floor(Int64,iss) round(distFromBlast,digits=2)   round(patient.SimParams[1],digits=2)  round(patient.SimParams[2],digits=2)    patient.SimParams[3] round(patient.SimParams[4],digits=2) patient.mobile bleed xLat xLon hitLat hitLon])
               end
            end
            sscore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
           # @yield request(patient.transport) #resource to let the victim function continue after transport
            @process victimEvac(env, patient, outputDB, mediclist, boatList, ambulist,medevac, transportqueue,transportqueue1, paramset,treatment_data,fmpUTM,homeUTM,ccpUTM,fmp_LLA,ccp_LLA,timeCCPHome,timeccpFMP,timeFMPHome,consts,inputDB,m)

        end
        if (event.cbrn == 4) #nuclear generation
            for j in 1:consts.numberofVictims
                    
                    
                x0 = victimPos[i,1]
      
       
                y0 = victimPos[i,2]
       
                hitloc = (victimPos[j,3],victimPos[j,4])
                
                distFromBlast = sqrt((hitloc[1]-x0)^2+(hitloc[2]-y0)^2)
                    
                    
            end
               
                iss = 0
                ip = 7 
               
                triage = 2
                #triage = getTriage(iss,ip)
                if iss > 0
                    timeOfDeath = 43500*iss^(-1.95)
                else 
                    timeOfDeath = 1.0e8
                   
                end
              
                SimParams=[]
                treatmentParams=[0.0,0.0,0.0,0.0,0.0]
                
                hospitalNeed = injury.hospitalNeed
                armor = 100.0
                scheduledTime = 0.0
                facility="PoI"
                patient = Patient(env, victag, age, is_treated, canwalk , lethally_injured,triage,x0,y0,iss,timeOfDeath,scheduledTime,SimParams,maxSimScore,bleed,(0,0),injury,injury2,injury3,hospitalNeed,consts,SS_array,distFromBlast,treatmentParams,threat,facility,ip)
                patient.hospitalneed=updateHospitalNeed(patient)
                if (lethally_injured == 1)
                    patient.SimParams[3]=1.0
                end
                
                
               realToD = 20000.0

                for i in 1:2000
                    test = checkHealth(env,patient,convert(Float64,i),0.01)
                    if test == 1
                        realToD = convert(Float64,i)
                        break
                    end
                end
                       
                if consts.creator
                 #   push!(victims,patient)
                  if threat == 0
                    push!(victims_list0, [victag floor(Int64,age) triage event.name round(hitloc[1],digits=2) round(hitloc[2],digits=2)  round(patient.x0,digits=2)  round(patient.y0,digits=2) round(distFromBlast,digits=2) floor(Int64,iss)  round(patient.SimParams[1],digits=2)  round(patient.SimParams[2],digits=2)    patient.SimParams[3] round(patient.SimParams[4],digits=2) patient.mobile bleed round(patient.SimParams[6],digits=2) round(patient.SimParams[7],digits=2) round(patient.SimParams[8],digits=2) round(patient.SimParams[9],digits=2) patient.dose])
                  end
                  if threat == 1
                    push!(victims_list1, [victag floor(Int64,age) triage event.name round(hitloc[1],digits=2) round(hitloc[2],digits=2)  round(patient.x0,digits=2)  round(patient.y0,digits=2) round(distFromBlast,digits=2) floor(Int64,iss)  round(patient.SimParams[1],digits=2)  round(patient.SimParams[2],digits=2)    patient.SimParams[3] round(patient.SimParams[4],digits=2) patient.mobile bleed round(patient.SimParams[6],digits=2) round(patient.SimParams[7],digits=2) round(patient.SimParams[8],digits=2) round(patient.SimParams[9],digits=2) patient.dose])
                  end
                  if threat == 2
                    push!(victims_list2, [victag floor(Int64,age) triage event.name round(hitloc[1],digits=2) round(hitloc[2],digits=2)  round(patient.x0,digits=2)  round(patient.y0,digits=2) round(distFromBlast,digits=2) floor(Int64,iss)  round(patient.SimParams[1],digits=2)  round(patient.SimParams[2],digits=2)    patient.SimParams[3] round(patient.SimParams[4],digits=2) patient.mobile bleed round(patient.SimParams[6],digits=2) round(patient.SimParams[7],digits=2) round(patient.SimParams[8],digits=2) round(patient.SimParams[9],digits=2) patient.dose])
                  end
                end
                sscore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
                @yield request(patient.transport) #resource to let the victim function continue after transport
                @process victim(env, patient, outputDB, firefighters,unit,amst, mediclist, ambulist,medevac, transportqueue, transportqueueT3, hospitalQueue, paramset,treatment_data,fmpUTM,homeUTM,ccpUTM,fmp_LLA,ccp_LLA,timeCCPHome,timeccpFMP,timeFMPHome,consts,inputDB,m,offroad)

            
        end
        
        
       if sscore <= 0.01 && patient.triage != 9

               patient.triage = 9
               VictimLogOutputDB(outputDB, patient.vicid,  now(env), 188,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,patient.hospitalID,patient.SimParams[3],patient.bloodvolume)
        else
            if consts.logg
                VictimLogOutputDB(outputDB, patient.vicid,  now(env), 42,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,patient.hospitalID,patient.SimParams[3],patient.bloodvolume)
            end
       end
       
      

   end
   if consts.creator
        if threat == 0    
            writedlm("victims0.txt",victims_list0)
        end
        if threat == 1    
            writedlm("victims1.txt",victims_list1)
        end
        if threat == 2    
            writedlm("victims2.txt",victims_list2)
        end
        
    end
      @process lastVictim(env, firefighters, mediclist, ambulist,medevac, paramset.Policy,paramset,outputDB,consts)

end
function logData(env::Environment,patient::Patient,action::Int64,paramset,consts,outputDB)
   #not really logging more like probing state
   
        
       sscore=generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts) 
       sscore = round(sscore,digits=3)
       # standing=checkMobility(env,patient)

       # if standing == true
       #     intStanding = 1
       # else
       #     intStanding = 0
       # end

       
       if sscore <= 0.0001 && patient.isdead == 0
           
           patient.triage = 9
           patient.isdead = 1
        
        
           VictimLogOutputDB(outputDB, patient.vicid,  now(env), 188,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,patient.hospitalID,patient.SimParams[3],patient.bloodvolume)

       end
       if patient.triage != 9
         
           if consts.logg && action!=1004
               VictimLogOutputDB(outputDB, patient.vicid,  now(env), action,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,patient.hospitalID,patient.SimParams[3],patient.bloodvolume)
           end
           if consts.logg && action==1004
            
                traveltime = timerandomised(4000/(8.0*60), "Travel",consts) #4km/15.33 m/s time in min
             
                VictimLogOutputDB(outputDB, patient.vicid,  now(env)+traveltime , action,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,patient.hospitalID,patient.SimParams[3],patient.bloodvolume)
           end

           if consts.logg && action==1004
            
                traveltime = timerandomised(4000/(8.0*60), "Travel",consts) #4km/15.33 m/s time in min
             
                VictimLogOutputDB(outputDB, patient.vicid,  now(env)+traveltime , action,sscore,patient.x0,patient.y0,patient.triage,patient.mobile,patient.isBleeding,paramset,patient.hospitalID,patient.SimParams[3],patient.bloodvolume)
           end

       end
   
end
@resumable function smallNoria(env::Environment, patient ::Patient,ambulist::AmbuList,paramset,timeccpFMP,consts,outputDB,fmpUTM)
    if paramset.TriageLevel == "None"
        priorityScore = 10
    else
        priorityScore = vicPriorityScore(env,patient,consts)
    end
    @yield request(ambulist.ForwardMedevac; priority = priorityScore)
    ResourceLogOutputDB(outputDB, "ForwardMedevac",ambulist.ForwardMedevac.capacity, ambulist.ForwardMedevac.level, now(env), "Request small noria",paramset)
    if isdead(patient)
        #morgue
        @yield release(ambulist.ForwardMedevac)
        ResourceLogOutputDB(outputDB, "ForwardMedevac",ambulist.ForwardMedevac.capacity, ambulist.ForwardMedevac.level, now(env), "release dead",paramset)
        return
    end
    logData(env,patient,14,paramset,consts,outputDB)
    
   
    if timeccpFMP !=0
 
        @yield timeout(env, timerandomised(timeccpFMP,"Other",consts))
        patient.x0 = fmpUTM.x
        patient.y0 = fmpUTM.y
        logData(env,patient,15,paramset,consts,outputDB)
      

        #@yield timeout(env, timerandomised(TotalSmallNoriaTransportTime,"Other"))
       @yield release(ambulist.ForwardMedevac)
       ResourceLogOutputDB(outputDB, "ForwardMedevac",ambulist.ForwardMedevac.capacity, ambulist.ForwardMedevac.level, now(env), "release small noria not offroad",paramset)
    # ambulist.TacMedevac.capacity += 1
    else #offroad
        patientUTM = Geodesy.UTM(patient.x0,patient.y0,0.0)
        distToFMP = getDistance(patientUTM,fmpUTM) #vol d'oiseau

        timeccpFMP = distToFMP/(30.0*60.0)+consts.UnloadFMP#30km/h suppose loading = unloading time
        @yield timeout(env, timerandomised(timeccpFMP,"Other",consts))
        patient.x0 = fmpUTM.x
        patient.y0 = fmpUTM.y
        logData(env,patient,15,paramset,consts,outputDB)
        @yield release(ambulist.ForwardMedevac)
        ResourceLogOutputDB(outputDB, "ForwardMedevac",ambulist.ForwardMedevac.capacity, ambulist.ForwardMedevac.level, now(env), "release small noria",paramset)

    end
  
end
@resumable function victim(env::Environment, patient ::Patient, outputDB ::SQLite.DB, firefighters ::Resource,unit::DeconUnit,amst::Amst, mediclist ::MedicList, ambulist ::AmbuList,medevac::Medevac, transportqueue ::Array, transportqueueT3 ::Array, hospitalQueue ::Array, paramset ::Parameterset,treatment_data::Array,fmpUTM::UTM,homeUTM::UTM,ccpUTM::UTM,fmp_LLA,ccp_LLA,timeCCPHome,timeccpFMP,timeFMPHome,consts,inputDB,m,offroad)
    
    # when condition is walking
    if isdead(patient)
        return
    end
 
    if patient.mobile == 1
        if paramset.Tourniquet
            if patient.isBleeding == 1
                tourniquet_success = rand()
                if tourniquet_success >= 0.3 #70% success rate of TQ application
                   
                    patient.isBleeding = 0
                    patient.SimParams[3] = 0.97 #not bleeding will survive
                    @yield request(patient.tourniquet)
                    @yield timeout(env, timerandomised(1.0, "Rescue",consts)) #♦time to apply tourniquet
                    logData(env,patient,44,paramset,consts,outputDB)
                    
                else
                    @yield timeout(env, timerandomised(1.0, "Rescue",consts)) #♦time to apply tourniquet
                    println("tq unsuccessful for $(patient.vicid)")
                end
                
                
            end
        end
        
        if patient.SimParams[5] < 7 && patient.SimParams[5] > 3
            #apply combopen
               combopen_success = rand()
               if combopen_success >= 0.1 #70% success rate of combopen application
                 if patient.isBleeding == 0
                   patient.SimParams[3] = 0.97 #if not bleeding will survive
                 end
                patient.SimParams[5] -= 1
                @yield request(patient.combopen)
                @yield timeout(env, timerandomised(0.2, "Rescue",consts)) #♦time to apply combopen and for helper to reach you
                logData(env,patient,31,paramset,consts,outputDB)
                                      
               else
                println("combopen unsuccessful for $(patient.vicid)")
               end
        
        end
        if patient.facility == "PoI" #only if spawning at PoI
            #self evacuate
        
            patientUTM = Geodesy.UTM(patient.x0,patient.y0,0.0)
            distToCCP = getDistance(patientUTM,ccpUTM) #vol d'oiseau

            timeToCCP = distToCCP/(consts.walkSpeed*60.0) #(5km/h) then converted to mins!
            
            @yield timeout(env, timerandomised(timeToCCP, "Rescue",consts))
            
            #triage
            if paramset.TriageLevel == "Nato"
                if patient.triage == 5
                
                    
                    
                    @yield timeout(env, timerandomised(timeCCPHome, "Rescue",consts))
                    patient.x0 = homeUTM.x
                    patient.y0 = homeUTM.y
                    logData(env,patient,0,paramset,consts,outputDB)
                
                    return
                
                else
                    if isdead(patient)
                        return
                    end
                                
                    patient.x0 = ccpUTM.x
                    patient.y0 = ccpUTM.y
                    
                    logData(env,patient,1,paramset,consts,outputDB)
                    standingTriage = @process triageAtLocation(env, patient, mediclist,"CCP",consts.TriageTimeUrgentDefault,paramset,consts,outputDB)
                    @yield standingTriage
          
                    
                end
            end          
        end
    else #immobile
        if paramset.Tourniquet
            if patient.isBleeding == 1 #apply tourniquet
         
         
                help,tag = scoutForMedic((patient.x0,patient.y0),5.0)
            
                if help == 1 && patient.tourniquet.capacity > 0
                    patient.SimParams[3] = 0.97
                    @yield request(patient.tourniquet)
                    patient.tourniquet.capacity -=1
                    @yield timeout(env, timerandomised(3.0, "Rescue",consts)) #♦time to apply tourniquet and for helper to reach you
                    logData(env,patient,44,paramset,consts,outputDB)
                    patient.isBleeding = 0
            
                else
        
                end
            end
           
           
        end
       
        if patient.SimParams[5] < 7 && patient.SimParams[5] > 3
            #apply combopen
               combopen_success = rand()
               if combopen_success >= 0.3 #70% success rate of combopen application
                if patient.isBleeding == 0
                   patient.SimParams[3] = 0.97 #not bleeding will survive
                end 
                patient.SimParams[5] -= 1
                @yield request(patient.combopen)
                patient.combopen.capacity -=1
                @yield timeout(env, timerandomised(0.2, "Rescue",consts)) #♦time to apply combopen and for helper to reach you
                logData(env,patient,31,paramset,consts,outputDB)

                   
               else
             
               end
        
        end

        if patient.facility ==  "PoI" && consts.inputDBname=="ukr.sqlite"
            # Search and Rescue by firefighters
            @yield timeout(env,rand()*0.01) 
            logData(env,patient,2,paramset,consts,outputDB) #waititng to be rescued
            
            # searchandrescuedirect = @process searchAndRescueUrgent(env, patient, firefighters, mediclist,consts,ccpUTM,paramset,outputDB)
            # @yield searchandrescuedirect
            logData(env,patient,18,paramset,consts,outputDB)
                # time-out for rescue duration
               
         
                patientUTM = Geodesy.UTM(patient.x0,patient.y0,0.0)
                distToCCP = getDistance(patientUTM,ccpUTM) #vol d'oiseau
        
               
                timeToCCP = distToCCP/(1.38*60.0) #(5km/h) then converted to mins!
                
                @yield timeout(env, timerandomised(timeToCCP, "Rescue",consts))
                patient.x0 = ccpUTM.x
                patient.y0 = ccpUTM.y
                logData(env,patient,1,paramset,consts,outputDB)
        # Triage
        # In case dead -> morgue, otherwise
            if isdead(patient)
                #morgue
                logData(env,patient,188,paramset,consts,outputDB)
                return
            else
                
                if paramset.TriageLevel == "Nato"
                        # triageprocess = @process triageAtLocation(env, patient, mediclist,"CCP",TriageTimeUrgentDefault,paramset)
                        # @yield triageprocess
                    if patient.triage != 9
            
                        patient.x0 = ccpUTM.x
                        patient.y0 = ccpUTM.y
                        
                        standingTriage = @process triageAtLocation(env, patient, mediclist,"CCP",consts.TriageTimeUrgentDefault,paramset,consts,outputDB)
                        @yield standingTriage
                        
                       
                    end
                end

            
            end
        end
        if patient.facility ==  "PoI" && consts.inputDBname!="ukr.sqlite"
            # Search and Rescue by firefighters
            @yield timeout(env,rand()*0.01) 
            logData(env,patient,2,paramset,consts,outputDB) #waititng to be rescued
            
            searchandrescuedirect = @process searchAndRescueUrgent(env, patient, firefighters, mediclist,consts,ccpUTM,paramset,outputDB)
            @yield searchandrescuedirect
        # Triage
        # In case dead -> morgue, otherwise
            if isdead(patient)
                #morgue
                logData(env,patient,188,paramset,consts,outputDB)
                return
            else
                
                if paramset.TriageLevel == "Nato"
                        # triageprocess = @process triageAtLocation(env, patient, mediclist,"CCP",TriageTimeUrgentDefault,paramset)
                        # @yield triageprocess
                    if patient.triage != 9
            
                        patient.x0 = ccpUTM.x
                        patient.y0 = ccpUTM.y
                        
                        standingTriage = @process triageAtLocation(env, patient, mediclist,"CCP",consts.TriageTimeUrgentDefault,paramset,consts,outputDB)
                        @yield standingTriage
                        
                       
                    end
                end

            
            end
        end
                
    end

    #...Stay and Play
    if paramset.Policy == "StayPlay"
        
        #...Treatment
        
        if isdead(patient)
            #morgue
            #ólogData(env,patient,3)
            return
        end
        if patient.contaminated != 0 && patient.facility == "PoI"
          
            fmp = Geodesy.LLA(47.23570323, 18.16913678) #R1 SVN/USA force transfer for decon!
            utm_location = UTMfromLLA(consts.quadrant,true,wgs84)
            fmpU = utm_location(fmp)
            smallnoria = @process smallNoria(env, patient,ambulist,paramset,timeccpFMP,consts,outputDB,fmpU)
            @yield smallnoria
            patient.x0 = fmpU.x
            patient.y0 = fmpU.y
            #start Decon
            
            @process decontamination(env, patient, unit,consts,paramset,outputDB)

        end
      
        if patient.triage == 5 # CCP -> NUCA
            # NUCA transport
              
               
                timeToCCP = timeCCPHome #(5km/h) then converted to mins!
                
                patient.x0 = homeUTM.x
                patient.y0 = homeUTM.y
                @yield timeout(env, timerandomised(timeToCCP, "Rescue",consts))
                logData(env,patient,7,paramset,consts,outputDB)
                
                
                return
            
            
        else
           
            if patient.facility == "PoI" && patient.contaminated == 0
            #logData(env,patient,14)
                smallnoria = @process smallNoria(env, patient,ambulist,paramset,timeccpFMP,consts,outputDB,fmpUTM)
                @yield smallnoria
            end
           
        end

        # Arrived FMP
       
        if isdead(patient)
            logData(env,patient,6,paramset,consts,outputDB)
            #morgue
            return
        else
        # treatment fmp
           
            treatmentprocess = @process treatmentAtLocation(env, patient, mediclist, paramset,"FMP",treatment_data,consts,outputDB)
            @yield treatmentprocess

        end
        #end
       
        #...Transport from FMP to Hosp

        # after treatment correct triage is known
        if isdead(patient)
            #morgue
            logData(env,patient,188,paramset,consts,outputDB)
            #logData(env,patient,43)
            return
        end

        if patient.triage == 5 || patient.triage ==3# FMP -> NUCA
            # home or NUCA transport
                timeFMPHome = 30.0 #irrelevant for mortality just for realism
                @yield timeout(env, timerandomised(timeFMPHome,"Other",consts))
                patient.x0 = homeUTM.x
                patient.y0 = homeUTM.y
                logData(env,patient,7,paramset,consts,outputDB)
                
                return
          
        else
        # triage is 1,2,4
           
          
            sscore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
            score = patient.triage*100+trunc(Int64,sscore)
            if paramset.TriageLevel == "None"
                score = 10
            end
            #@yield request(ambulist.transportcheckoccupied; priority = 1)
            
            enterVicTranspQueue!(transportqueue, patient, score)
            
            locationx=patient.x0
            locationy=patient.y0
            #callin to transport
            @yield timeout(env, timerandomised(2.0,"Other",consts)) #timeout to call in the ambulance
            logData(env,patient,16,paramset,consts,outputDB)
            # @yield release(ambulist.transportcheckoccupied)
            transp = @process transport(env, transportqueue, transportqueueT3, mediclist, ambulist, hospitalQueue, paramset,treatment_data,consts,inputDB,outputDB,locationx,locationy,m,offroad,homeUTM)                        
            

            @yield transp
            @yield request(patient.transport)  
           
        end
        
                        
    end #stayplay end
            
    #...Scoop and Run
    
    if paramset.Policy == "ScoopRun"
        #...Transport
      
        if isdead(patient)
            #morgue
           # logData(env,patient,3)
            return
        end
        
        
       
        if triageMistakes(patient.triage, patient.triagecorrectness) == 3 #T3 victims
              # CCP -> home or NUCA transport

            
            @yield timeout(env, timerandomised(timeCCPHome,"Other",consts))
            patient.x0 = homeUTM.x
            patient.y0 = homeUTM.y
       
            logData(env,patient,0,paramset,consts,outputDB)
                return
                                    
        end
        if triageMistakes(patient.triage, patient.triagecorrectness) in [0,1,2,4]
                
                # triaged as T1,2,4
                if paramset.TriageLevel == "None"
                    score = 10
                    
                else
                
                    sscore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
                  
                    score = triageMistakes(patient.triage, patient.triagecorrectness)*100+trunc(Int64,sscore)
                    
                end
                #@yield request(ambulist.transportcheckoccupied; priority = 1)


                enterVicTranspQueue!(transportqueue, patient, score)
                
                
                patient.x0 = ccpUTM.x
                patient.y0 = ccpUTM.y
                @yield timeout(env, timerandomised(2.0,"Other",consts)) #timeout to call in the ambulance
                logData(env,patient,16,paramset,consts,outputDB)

               # @yield release(ambulist.transportcheckoccupied)
                locationx=patient.x0
                locationy=patient.y0             
                transp1 = @process transport(env, transportqueue, transportqueueT3, mediclist, ambulist, hospitalQueue, paramset,treatment_data,consts,inputDB,outputDB,locationx,locationy,m,offroad,homeUTM)
                @yield transp1
                @yield request(patient.transport)   
                
        else #triage 9
          #  @yield release(ambulist.transportcheckoccupied)
            return
        end
    
    end
    
end #victim function
@resumable function dischargePatient(env ::Environment, patient::Patient,hospitalQueue ::Array{Hospital,1}, hospitalID ::Int64,paramset::Parameterset,consts,outputDB)
    for hospital in hospitalQueue
       
       
        if hospital.ID == hospitalID
            if patient.triage == 1 && patient.age > 15
                #hospital.CapT1 += 1
                
                #NlogData(env,patient,26,paramset,consts,outputDB)
                
                #NlogData(env,patient,26,paramset,consts,outputDB)
                patient.is_treated=1
                
                return
                
            end
            if patient.triage == 2 && patient.age > 15
                #hospital.CapT2 += 1
         
                #logData(env,patient,26,paramset,consts,outputDB)
                #logData(env,patient,26,paramset,consts,outputDB)
                patient.is_treated=1
               
                return
            end
            if patient.triage == 3 && patient.age > 15
           
               # hospital.CapT3 += 1
                #logData(env,patient,26,paramset,consts,outputDB)
                #logData(env,patient,26,paramset,consts,outputDB)
                patient.is_treated=1
                
                return
            end

            if patient.triage == 1 && patient.age <= 15
              #  hospital.CapT1Ped += 1
             
                #logData(env,patient,26,paramset,consts,outputDB)
                #logData(env,patient,26,paramset,consts,outputDB)
                patient.is_treated=1
              
                return
            end
            if patient.triage == 2 && patient.age <= 15
              #  hospital.CapT2Ped += 1
         
               # logData(env,patient,26,paramset,consts,outputDB)
               # logData(env,patient,26,paramset,consts,outputDB)
                patient.is_treated=1
                
                return
            end
            if patient.triage == 3 && patient.age <= 15
               # hospital.CapT3Ped += 1
         
               # logData(env,patient,26,paramset,consts,outputDB)
               # logData(env,patient,26,paramset,consts,outputDB)
                patient.is_treated=1
                return
            end
        end
    end
    
end
function checkDestinationHospital(patient::Patient,inputDB::SQLite.DB,consts::Constants)

    id = patient.hospitalID
    hospLLAs = DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Hospital" ))


    for i in 1:size(hospLLAs,1)
        if id == hospLLAs[i,:].ID
            hospX = hospLLAs[i,:].LAT
            hospY = hospLLAs[i,:].LON
            hosp_LLA = Geodesy.LLA(hospX,hospY) #HopsID UTM coords
            #•ospUTM = utm_location(hosp_LLA)
                
                utm_location = UTMfromLLA(consts.quadrant,true,wgs84)
                lla_location = LLAfromUTM(consts.quadrant,true,wgs84)
                originLLA = lla_location(Geodesy.UTM(patient.x0,patient.y0,0.0))
                hospUTM=utm_location(hosp_LLA)
                originECF = ECEFfromLLA(wgs84)(originLLA)
                hospECF = ECEFfromLLA(wgs84)(hosp_LLA)
            
                distToHosp=OpenStreetMapX.distance(originECF.x,originECF.y,originECF.z,hospECF.x,hospECF.y,hospECF.z)
                if distToHosp < 2.0
                    #no need to transport
                    return false
                else 
                    return true

                end
        end
    end

end
@resumable function searchAndRescueUrgent(env::Environment, patient :: Patient, firefighters :: Resource, mediclist :: MedicList,consts::Constants,ccpUTM::UTM,paramset::Parameterset,outputDB)

    
    #     @yield (requestfirefighters = request(firefighters; priority = 500)) 
    # else
    #     priorityScore= 500
    #     @yield (requestfirefighters = request(firefighters; priority = priorityScore ))
    # end

    # if state(requestfirefighters) == SimJulia.processed
        

        if paramset.PreTriage
        
           #2 @process releaseMedicDirect(env, firefighters, requestfirefighters)
            
            # time-out for pr-triage duration
            @yield timeout(env, timerandomised(consts.PreTriageTime,"Normal",consts))
            

            if isdead(patient)
                #morgue
                logData(env,patient,188,paramset,consts,outputDB)         
                return
            else

                priorityScore = vicPriorityScore(env,patient,consts)
                
                @yield request(firefighters; priority = priorityScore)
                ResourceLogOutputDB(outputDB, "Firefighters",firefighters.capacity, firefighters.level, now(env), "Request",paramset)
                logData(env,patient,18,paramset,consts,outputDB)
                # time-out for rescue duration
               
         
                patientUTM = Geodesy.UTM(patient.x0,patient.y0,0.0)
                distToCCP = getDistance(patientUTM,ccpUTM) #vol d'oiseau
        
               
                timeToCCP = distToCCP/(1.38*60.0) #(5km/h) then converted to mins!
                
                @yield timeout(env, timerandomised(timeToCCP, "Rescue",consts))
                patient.x0 = ccpUTM.x
                patient.y0 = ccpUTM.y
                logData(env,patient,1,paramset,consts,outputDB)
                @yield release(firefighters)
                
            end
        
        else #noPreTriage
            
            if isdead(patient)
                #morgue
                logData(env,patient,188,paramset,consts,outputDB)      
                return
            else
                
                
            # time-out for rescue duration
                
                priorityScore = randomSort()
                @yield request(firefighters; priority = priorityScore)
                logData(env,patient,18,paramset,consts,outputDB)
        
                patientUTM = Geodesy.UTM(patient.x0,patient.y0,0.0)
                distToCCP = getDistance(patientUTM,ccpUTM)
        
              
                timeToCCP = distToCCP/(1.38*60.0) #(5km/h) then converted to mins! Vol d'oiseau
                
                @yield timeout(env, timerandomised(timeToCCP, "Rescue",consts))
                patient.x0 = ccpUTM.x
                patient.y0 =ccpUTM.y
                logData(env,patient,1,paramset,consts,outputDB)
                @yield release(firefighters)
            end
        
        end
    # end
         
end
@resumable function triageAtLocation(env::Environment, patient :: Patient, mediclist :: MedicList,loc ::String,triageTime::Float64,paramset::Parameterset,consts,outputDB)

    if loc == "CCP"
        triageDoctor = mediclist.ccpDoctor
        triageNurse = mediclist.ccpNurse
                    
    elseif loc == "FMP"
        triageDoctor = mediclist.fmpDoctor
        triageNurse = mediclist.fmpNurse
      
              
    else

    end
    @yield (requesttriagedoctor = request(triageDoctor; priority = 2)) | (requesttriagenurse = request(triageNurse; priority = 2))
    if state(requesttriagedoctor) == SimJulia.processed
        triageMed = "d"
        if state(requesttriagenurse) == SimJulia.processed
            @yield release(mediclist.ccpNurse)
        else
            @process releaseMedicDirect(env, triageNurse, requesttriagenurse)
        end
    elseif state(requesttriagenurse) == SimJulia.processed
        triageMed = "n"
        @process releaseMedicDirect(env, triageDoctor, requesttriagedoctor)
    end
    
    # time-out for pr-triage duration
    @yield timeout(env, timerandomised(triageTime,"Other",consts))
    patient.triage = setTriage(env,patient,consts)
    @yield timeout(env, timerandomised(consts.TriageTimeUrgentDefault,"Other",consts)) #triage time
    logData(env,patient,19,paramset,consts,outputDB)
    # release correct medic
    if triageMed == "d"
        release(triageDoctor)
    elseif triageMed == "n"
        release(triageNurse)
    
    end
end
@resumable function transport(env ::Environment, transportqueue ::Array, transportqueueT3 ::Array, mediclist ::MedicList, ambulist ::AmbuList, hospitalQueue ::Array{Hospital,1}, paramset ::Parameterset,treatment_data::Array,consts::Constants,inputDB,outputDB,locationx,locationy,m,offroad,homeUTM)
    
    if length(transportqueue) == 0 
      #  @yield release(ambulist.transportcheckoccupied)
  
  
        return
    end
          
    if length(transportqueue) != 0
        if isdead(transportqueue[1][2])
            deleteat!(transportqueue,1)
   
            return
        end
    end
   
    if (ambulist.TacMedevac.level > ambulist.TacMedevac.capacity) 
      
       # @yield release(ambulist.transportcheckoccupied)
      
        println("no more ambulances")
        return
           
    end

    # ambulances are present
    
    # dravailable = (mediclist.fmpMMT.capacity+mediclist.fmpDoctor.capacity-mediclist.nrDoctorTransport > mediclist.minDoctorFMP) && (mediclist.fmpMMT.level < mediclist.fmpMMT.capacity || mediclist.fmpDoctor.level < mediclist.fmpDoctor.capacity)
    # nuavailable = (mediclist.fmpMMT.capacity+mediclist.fmpNurse.capacity-mediclist.nrNurseTransport > mediclist.minNurseFMP) && (mediclist.fmpMMT.level < mediclist.fmpMMT.capacity || mediclist.fmpNurse.level < mediclist.fmpNurse.capacity)
    
    if mediclist.nrNurseTransport>0
        nuavailable = true
    else
        nuavailable = false
    end
    if mediclist.nrDoctorTransport>0
        dravailable = true
    else
        dravailable = false
    end


    
    if (((patient, supervision, hospid) = findVictimTransport(env,transportqueue, hospitalQueue, paramset.supervisionlevels, (dravailable,nuavailable), paramset.HospitalDistribution,inputDB,paramset,locationx,locationy,m,consts.quadrant,consts.offroad)) == (nothing,nothing,nothing))
        
        return
    end
    
    
    # All checks satisfied
    if  patient.triage != 9 && patient.triage != 3 && checkDestinationHospital(patient,inputDB,consts)
            @yield request(ambulist.ForwardMedevac)
            
          
           # ambulist.ForwardMedevac.level += 1 
            ResourceLogOutputDB(outputDB, "ForwardMedevac",ambulist.ForwardMedevac.capacity, ambulist.ForwardMedevac.level, now(env), "Request",paramset)
            
            
            @yield timeout(env, timerandomised(consts.LoadFMP,"Other",consts)) #timeout to call in the ambulance
            logData(env,patient,20,paramset,consts,outputDB)
        
    end 
    if !checkDestinationHospital(patient,inputDB,consts)
       # println("patient $(patient.vicid) is already at a R2 no transfer needed he needed $(patient.hospitalneed)")
    end
    if patient.triage == 3
        
        return
   

    end
    
    #~patient.transporthasstarted = true

    if paramset.Policy == "StayPlay" 
        # correct triage after treatment
        if patient.triage != 9
            triagecat = patient.triage
            # @yield timeout(env, timerandomised(consts.TriageTimeUrgentDefault,"Other",consts)) #triage time
            # logData(env,patient,19,paramset,consts,outputDB)

        else
            triagecat = 9
        end
    elseif paramset.Policy == "ScoopRun"
        # incorrect triage
        if patient.triage != 9
            triagecat = triageMistakes(patient.triage,patient.triagecorrectness)
        else
            triagecat = 9
        end
        
    end
       
   
    if patient.triage != 9 
        requestHospital(env, hospitalQueue, hospid, triagecat, patient.age >= 15) #not actual request as hospital is not modelled as resource
        if hospid !== nothing
            patient.hospitalID=hospid
          
        else
            @yield timeout(env, 60.0)
            requestHospital(env, hospitalQueue, hospid, triagecat, patient.age >= 15) #not actual request as hospital is not modelled as resource
            if hospid !== nothing
                patient.hospitalID=hospid
              
            end
        end
    else
       
    
    end
      

    if  patient.triage != 9
       # @yield timeout(env,rand()*0.0001)
       
       transportScore = generateSimedisScoreGompertzPatient(env,(now(env)-patient.scheduledTime)*60,patient,consts)
        if paramset.mascal
            if transportScore < 5.0
                patient.triage = 4
                logData(env,patient,21,paramset,consts,outputDB)
                
              #  @yield release(ambulist.transportcheckoccupied)
                @yield release(ambulist.TacMedevac) 
                return
            end
                logData(env,patient,22,paramset,consts,outputDB)
      
        end
    end
 

    loadtime = timerandomised(consts.LoadFMP, "Other",consts)
    traveltime = timerandomised(findHospTravelTime(hospitalQueue, hospid,patient,m,consts.quadrant), "Travel",consts)
   
    
       #medicskill = if supervision == 0 1 else supervision end


   # EMT treatment

   treatStartTime = trunc(Int64,now(env))
    attack=rand()
    if paramset.Policy == "StayPlay" 
        hit = consts.ambushRate/6.0
    end
    if paramset.Policy == "ScoopRun" 
        hit = consts.ambushRate
    end
    @yield timeout(env, timerandomised(consts.LoadFMP,"Other",consts)) #load patient to ambulance
    logData(env,patient,23,paramset,consts,outputDB)
    patient.treatmentParams[2]=0.5
    patient.treatmentParams[3]=20.0
    patient.treatmentParams[4]= (now(env)-patient.scheduledTime)*60
    patient.treatmentParams[5]=(((10.0+patient.iss)/patient.triage)+(now(env)-patient.scheduledTime)*60)
      
    if attack >= hit
        
        @yield timeout(env, (loadtime + traveltime + consts.UnloadHosp)/2.0)
        if patient.is_treated == 0 && patient.triage != 9
            logData(env,patient,24,paramset,consts,outputDB)
        end
        @yield timeout(env, (loadtime + traveltime + consts.UnloadHosp)/2.0)
        logData(env,patient,25,paramset,consts,outputDB)
        
        #release(patient.transport)
    
        @yield timeout(env, consts.HospDropOff)
        @yield release(ambulist.ForwardMedevac)
        ResourceLogOutputDB(outputDB, "ForwardMedevac",ambulist.ForwardMedevac.capacity, ambulist.ForwardMedevac.level, now(env), "Released",paramset)
        #put arrival hospital inside transport end
        mapHospital = patient.hospitalID #location of the patient
        hospLLAs = DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Hospital" ))

        for i in 1:size(hospLLAs,1)
            if mapHospital == hospLLAs[i,:].ID

                if !offroad
                    hospX = hospLLAs[i,:].LAT
                    hospY = hospLLAs[i,:].LON
                    hosp_LLA = Geodesy.LLA(hospX,hospY) #HopsID UTM coords
                    utm_location = UTMfromLLA(consts.quadrant,true,wgs84)
                    hospUTM = utm_location(hosp_LLA)
                    hospNode = point_to_nodes((hospX,hospY),m)
                    fmpNode = point_to_nodes((fmp_LLA.lat,fmp_LLA.lon),m)
                    route, distance,route_time = OpenStreetMapX.fastest_route(m,hospNode,fmpNode)
                    timeToFMP = (route_time/60.0)*1.7
                    @yield timeout(env, timerandomised(timeToFMP,"Other",consts))
                    patient.x0 = hospUTM.x
                    patient.y0 = hospUTM.y
                else
                    hospX = hospLLAs[i,:].LAT
                    hospY = hospLLAs[i,:].LON
                    hosp_LLA = Geodesy.LLA(hospX,hospY) #HopsID UTM coords
                    #•ospUTM = utm_location(hosp_LLA)
                    
                    utm_location = UTMfromLLA(consts.quadrant,true,wgs84)
                    lla_location = LLAfromUTM(consts.quadrant,true,wgs84)
                    originLLA = lla_location(Geodesy.UTM(patient.x0,patient.y0,0.0))
                    hospUTM=utm_location(hosp_LLA)
                    originECF = ECEFfromLLA(wgs84)(originLLA)
                    hospECF = ECEFfromLLA(wgs84)(hosp_LLA)
                
                    distToHosp=OpenStreetMapX.distance(originECF.x,originECF.y,originECF.z,hospECF.x,hospECF.y,hospECF.z)
                    driveTime = (distToHosp*2.0/(30*(16.6666))) 
                    
                    timeToHosp = (driveTime/60.0)*1.7
                    
                    @yield timeout(env, timerandomised(timeToHosp,"Other",consts))
                    patient.x0 = hospUTM.x
                    patient.y0 = hospUTM.y
                                    
                end
                                        
            else
                
            end
        end
            
        logData(env,patient,17,paramset,consts,outputDB) 
        if length(hospitalQueue)>0
                                                        
            # add a patient.hospital container for direction
            treatmenthospurgent=@process treatmentInHospital(env, patient,hospitalQueue,patient.hospitalID,treatment_data,paramset,consts,outputDB,homeUTM)
            @yield treatmenthospurgent
        end
    else
        logData(env,patient,187,paramset,consts,outputDB)  #ambush kills patient
       
        logData(env,patient,188,paramset,consts,outputDB)  #ambush kills patient
       

    end

   
    # if paramset.Policy == "StayPlay"
    #     if supervision == 0
    #         release(mediclist.fmpMMT)
    #         ResourceLogOutputDB(outputDB, "fmpDoctor", mediclist.fmpMMT.capacity, mediclist.fmpMMT.level, now(env), "Return from transport",paramset)
    #         mediclist.nrDoctorTransport -= 1
    #         mediclist.nrNurseTransport -= 1
    #     elseif supervision == 1
    #         mediclist.nrDoctorTransport -= 1
    #         if mediclist.fmpNurse.level < mediclist.fmpNurse.capacity #free nurse
    #             mediclist.fmpNurse.capacity -= 1 # take one avaialble nurse to form MMT
    #             mediclist.fmpDoctor.capacity -= 1
    #             mediclist.fmpDoctor.level -= 1
    #             @process addMedic(env, mediclist.fmpMMT, 1)
    #             ResourceLogOutputDB(outputDB, "fmpMMT", mediclist.fmpMMT.capacity, mediclist.fmpMMT.level, now(env), "Form MMT",paramset)
    #             ResourceLogOutputDB(outputDB, "fmpDoctor", mediclist.fmpDoctor.capacity, mediclist.fmpDoctor.level, now(env), "Form MMT",paramset)
    #             ResourceLogOutputDB(outputDB, "fmpNurse", mediclist.fmpNurse.capacity, mediclist.fmpNurse.level, now(env), "Form MMT",paramset)
    #         else
    #             release(mediclist.fmpDoctor)
    #         end
    #     elseif supervision == 2
            
    #         mediclist.nrNurseTransport -= 1
    #         if mediclist.fmpDoctor.level < mediclist.fmpDoctor.capacity #free doctor
    #             mediclist.fmpDoctor.capacity -= 1 # take one avaialable doctor to form MMT
    #             mediclist.fmpNurse.capacity -= 1
    #             mediclist.fmpNurse.level -= 1
    #             @process addMedic(env, mediclist.fmpMMT, 1)
    #             ResourceLogOutputDB(outputDB, "fmpMMT", mediclist.fmpMMT.capacity, mediclist.fmpMMT.level, now(env), "Form MMT",paramset)
    #             ResourceLogOutputDB(outputDB, "fmpDoctor", mediclist.fmpDoctor.capacity, mediclist.fmpDoctor.level, now(env), "Form MMT",paramset)
    #             ResourceLogOutputDB(outputDB, "fmpNurse", mediclist.fmpNurse.capacity, mediclist.fmpNurse.level, now(env), "Form MMT",paramset)
    #         else
    #             release(mediclist.fmpNurse)
    #         end
            
    #     end
    # end

    # @yield release(ambulist.TacMedevac)
   
end
@resumable function truckTransport(env ::Environment, transportqueue ::Array, trucklist ::TruckList, paramset ::Parameterset,consts::Constants,outputDB,isReady::Bool,pickUp::UTM,dropOff::UTM)
    
    if isReady == true
            if length(transportqueue) == 0 
            #  @yield release(ambulist.transportcheckoccupied)
        
                return
            end
         
        
            if (boatlist.ForwardMedevac.level > boatlist.ForwardMedevac.capacity) 
            
            # @yield release(ambulist.transportcheckoccupied)
              
                return
                        
            end
            
            @yield request(boatlist.ForwardMedevac)
        
        for i in 1:min(length(transportqueue),consts.TruckCapacity) #BoatCapacity spots on 1 boat
            
           
            # else
                if i < length(transportqueue) && length(transportqueue) > 0
                    
                    patient = transportqueue[1][2]
                    
                    # if  patient.triage != 9
                    #     # @yield timeout(env,rand()*0.0001)
                    #     transportScore = generateSimedisScoreGompertzPatient(env,now(env)*60,patient)
                    #     if paramset.mascal
                    #         if transportScore < 5.0
                    #             patient.triage = 4
                    #             logData(env,patient,21,paramset,consts,outputDB)
                                
                    #         #  @yield release(ambulist.transportcheckoccupied)
                    #            # @yield release(ambulist.TacMedevac) 
                    #             return
                    #         end
                            
                    #         logData(env,patient,22,paramset,consts,outputDB)
                    
                    #     end
                    # end
    
                    dropOffPointLLA = Geodesy.LLA(-4.272567211685849, 15.293562984507624,0.0)
                
                    utmlocation = UTMfromLLA(33,true,wgs84)
                    dropPointUTM=utmlocation(dropOffPointLLA)
                            
                    patient.x0 = dropPointUTM.x
                    patient.y0 = dropPointUTM.y
                    logData(env,patient,23,paramset,consts,outputDB)
                    loadtime = timerandomised(0.2, "Other",consts)
                    traveltime = timerandomised(4000/(8.0*60), "Travel",consts) #4km/15.33 m/s time in min
                    logData(env,patient,1001,paramset,consts,outputDB)
                    @yield timeout(env, (loadtime + traveltime )) #UnloadBoat = 0.2 min
                    if length(transportqueue)>0
                        transportqueue = deleteat!(transportqueue,1) #take the 5 first patients
                    end
                   
                    release(patient.transport)
                
                end
                
            
            # end

        end

            boatlist.ForwardMedevac.capacity = boatlist.ForwardMedevac.capacity - 1
            
            @yield timeout(env, 0.3) #unloadFrom Boat
        
            @yield release(boatlist.ForwardMedevac)
            @yield timeout(env, timerandomised(4000/(8.0*60), "Travel",consts)) #return time
            boatlist.ForwardMedevac.capacity = boatlist.ForwardMedevac.capacity + 1
            
        else
          
    end
   
end
@resumable function boatTransport(env ::Environment, transportqueueX ::Array, boatlist ::BoatList, paramset ::Parameterset,consts::Constants,inputDB,outputDB,isReady::Bool,line::Int64)
    
    if isReady == true
            if length(transportqueueX) == 0
            #  @yield release(ambulist.transportcheckoccupied)
                
                return
            end    
        
            if (boatlist.ForwardMedevac.level > boatlist.ForwardMedevac.capacity) 
            
            # @yield release(ambulist.transportcheckoccupied)
                println("no more boats capa:$(boatlist.ForwardMedevac.capacity), level:$(boatlist.ForwardMedevac.level)")
                return
        
                
            end
            
            @yield request(boatlist.ForwardMedevac)
            dropOffPointLLA = Geodesy.LLA(-4.271781341476751, 15.293627469828829,0.0) 
            utmlocation = UTMfromLLA(33,true,wgs84)
            dropPointUTM=utmlocation(dropOffPointLLA)
       
        for i in 1:min(length(transportqueueX),consts.BoatCapacity) #BoatCapacity spots on 1 boat
               
                if length(transportqueueX) > 0
                  
                    
                    patient = transportqueueX[1][2] #♥patient is the first element (FIFO)
                    
                    
                    #push!(boatlist.patients,patient)
                    loadpatient = timerandomised(0.1, "Other",consts) #time to load 1 patient
                    @yield timeout(env, loadpatient)
                
                    # if  patient.triage != 9
                    #     # @yield timeout(env,rand()*0.0001)
                    #     transportScore = generateSimedisScoreGompertzPatient(env,now(env)*60,patient)
                    #     if paramset.mascal
                    #         if transportScore < 5.0
                    #             patient.triage = 4
                    #             logData(env,patient,21,paramset,consts,outputDB)
                                
                    #         #  @yield release(ambulist.transportcheckoccupied)
                    #            # @yield release(ambulist.TacMedevac) 
                    #             return
                    #         end
                            
                    #         logData(env,patient,22,paramset,consts,outputDB)
                    
                    #     end
                    # end
                    logData(env,patient,1001,paramset,consts,outputDB) #has been loaeded on boat
                    # for patient in boatlist.patients
                    #     patient.x0 = dropPointUTM.x
                    #     patient.y0 = dropPointUTM.y
                    # end
                    logData(env,patient,23,paramset,consts,outputDB)
                    patient.x0 = dropPointUTM.x
                    patient.y0 = dropPointUTM.y
                    #
                    if i == 1
                        loadtime = timerandomised(2.0, "Other",consts) #time to depart
                        traveltime = timerandomised(4000/(8.0*60), "Travel",consts) #4km/15.33 m/s time in min
                        @yield timeout(env, (loadtime + traveltime ))    
                    end
                       
                    logData(env,patient,1002,paramset,consts,outputDB)
                            
                    if length(transportqueueX) > 0
                        
                        deleteat!(transportqueueX,1) #take the BoatCapacity first patients
                      
                        release(patient.transport)
                    end
            
                        
                        
                else 
                            println("empty queue")  
                            return
                                            
                end
                   
                
            
         end

        ResourceLogOutputDB(outputDB, "zodiac", boatlist.ForwardMedevac.capacity, boatlist.ForwardMedevac.level, now(env), "dropOff",paramset)
        
        boatlist.ForwardMedevac.capacity = boatlist.ForwardMedevac.capacity - 1
        @yield release(boatlist.ForwardMedevac)
        @yield timeout(env, 0.3) #unloadFrom Boat
        
        #  
        @yield timeout(env, timerandomised(4000/(8.0*60), "Travel",consts)) #return time
        ResourceLogOutputDB(outputDB, "zodiac", boatlist.ForwardMedevac.capacity, boatlist.ForwardMedevac.level, now(env), "returned",paramset)
        # if length(transportqueueX) >= 0
        #     logData(env,patient,1004,paramset,consts,outputDB) #boat returned emplty from brazza
        # end
        boatlist.ForwardMedevac.capacity = boatlist.ForwardMedevac.capacity + 1
        return
         
    
          
    end
   
end
@resumable function treatmentInHospital(env ::Environment, patient ::Patient,hospitalQueue ::Array{Hospital,1}, hospitalID ::Int64,treatment_data::Array,paramset::Parameterset,consts,outputDB,homeUTM)

    patient.is_treated = 1 #set to 1 because the hospital treatment is the final step 
    
    logData(env,patient,29,paramset,consts,outputDB) #start hospital treatment
    patient.SimParams[3]=0.2
    applyTreatment(env,patient,((patient.iss/patient.triage)+5.0+now(env))*60,((patient.iss/patient.triage)+5.0)*60,3,now(env),treatment_data,paramset,consts)

    @yield timeout(env,(60.0))
    patient.x0 = homeUTM.x #sending patients to base
    patient.y0 = homeUTM.y
    logData(env,patient,30,paramset,consts,outputDB) #finish hosp treatment
    @process dischargePatient(env,patient,hospitalQueue,hospitalID,paramset,consts,outputDB)
    #@yield discharge
    
end
@resumable function treatmentAtLocation(env::Environment, patient ::Patient, mediclist ::MedicList, paramset ::Parameterset,loc::String,treatment_data::Array,consts,outputDB,transportqueue ::Array = Array{Tuple}(undef, 1))
    if paramset.TriageLevel == "None"
        priorityScore = 10
    else
        priorityScore = vicPriorityScore(env,patient,consts)
    end
    
    treatMed = ""
    medtype= ""
    
    # claim resources
    if loc == "CCP"
        treatmentDoctor = mediclist.ccpDoctor
        treatmentNurse = mediclist.ccpNurse
        elseif loc == "FMP"
            treatmentDoctor = mediclist.fmpDoctor
            treatmentNurse = mediclist.fmpNurse
   
        elseif loc == "On Site"
            treatmentDoctor = mediclist.onSiteDoctor
            treatmentNurse = mediclist.onSiteNurse
    else
            println("You need a location for Treatment!")
        
    end #check if location and adapt resources
    triagecat =setTriage(env,patient,consts)
    @yield timeout(env, timerandomised(consts.TriageTimeUrgentDefault,"Other",consts)) #triage time
    logData(env,patient,19,paramset,consts,outputDB) #start FMP treatment
    patient.triage = triagecat #patient is retriaged at FMP arrival
    # check if T1,T2,T4 victim is still alive
    if isdead(patient)
        @yield timeout(env, timerandomised(consts.TriageTimeUrgentDefault,"Other",consts))
        logData(env,patient,36,paramset,consts,outputDB) #found dead in treatment
        return
    end   
  
    
        patient.treatmenthasstarted = true
        
                    
        # time out for treatment duration
        if medtype == "fmpmmt"
            medicskill = 0
        elseif medtype == "fmpdr"
            medicskill = 1
        elseif medtype == "fmpnu"
            medicskill = 2
        end
       
        
        if loc == "FMP"
            
            @yield timeout(env, timerandomised(1.0,"Other",consts)) #handoff time to FMP
            logData(env,patient,38,paramset,consts,outputDB) #start FMP treatment
            patient.treatmentParams[1]=20.0
            patient.treatmentParams[2]=0.0
            patient.treatmentParams[3]=0.0
            patient.treatmentParams[4]= (now(env)-patient.scheduledTime)*60
            patient.treatmentParams[5]=((30.0)+(now(env)-patient.scheduledTime))*60
            if patient.isBleeding == 1 #fmp applies a tourniquet if needed !
                patient.SimParams[3]=0.9
            end
            treatmentTime = 30.0
        
            @yield timeout(env, treatmentTime )
         
            logData(env,patient,39,paramset,consts,outputDB) # FMP treatment ended
            triagecat =setTriage(env,patient,consts)
            @yield timeout(env, timerandomised(consts.TriageTimeUrgentDefault,"Other",consts)) #triage time
           
            @yield timeout(env, 1440.0 ) #timeout for evac in the strategyA
           
            logData(env,patient,19,paramset,consts,outputDB) #start FMP treatment
            if triagecat == 3
                logData(env,patient,7,paramset,consts,outputDB) #start FMP treatment
                return
            end
                
            
            
        end
          
     #finish urgent treatment

    # release resources (try to form mmt if possible)
    #= if medtype == "fmpmmt"
        @yield release(mediclist.fmpMMT)
    elseif  medtype == "fmpdr"
        if mediclist.fmpNurse.level < mediclist.fmpNurse.capacity #free nurse
           # println("taking a nurse to form an MMT...")
            mediclist.fmpNurse.capacity -= 1 # take one avaialble nurse to form MMT
            mediclist.fmpDoctor.capacity -= 1
            mediclist.fmpDoctor.level -= 1
            @process addMedic(env, mediclist.fmpMMT, 1)
            ResourceLogOutputDB(outputDB, "fmpMMT", mediclist.fmpMMT.capacity, mediclist.fmpMMT.level, now(env), "Form MMT",paramset)
            ResourceLogOutputDB(outputDB, "fmpDoctor", mediclist.fmpDoctor.capacity, mediclist.fmpDoctor.level, now(env), "Form MMT",paramset)
            ResourceLogOutputDB(outputDB, "fmpNurse", mediclist.fmpNurse.capacity, mediclist.fmpNurse.level, now(env), "Form MMT",paramset)
        else
            release(mediclist.fmpDoctor)
        end
    elseif medtype == "fmpnu"
        if mediclist.fmpDoctor.level < mediclist.fmpDoctor.capacity #free doctor
         #   println("taking a doctor to form MMT...")
            mediclist.fmpDoctor.capacity -= 1 # take one avaialable doctor to form MMT
            mediclist.fmpNurse.capacity -= 1
            mediclist.fmpNurse.level -= 1
            @process addMedic(env, mediclist.fmpMMT, 1)
            ResourceLogOutputDB(outputDB, "fmpMMT", mediclist.fmpMMT.capacity, mediclist.fmpMMT.level, now(env), "Form MMT",paramset)
            ResourceLogOutputDB(outputDB, "fmpDoctor", mediclist.fmpDoctor.capacity, mediclist.fmpDoctor.level, now(env), "Form MMT",paramset)
            ResourceLogOutputDB(outputDB, "fmpNurse", mediclist.fmpNurse.capacity, mediclist.fmpNurse.level, now(env), "Form MMT",paramset)
        else
            release(mediclist.fmpNurse)
        end
    end
    logData(env,patient,40)   =#
    
end
@resumable function lastVictim(env ::Environment, firefighters ::Resource, mediclist ::MedicList, ambulist ::AmbuList,medevac::Medevac, policy ::String,paramset::Parameterset,outputDB,consts)

    @yield timeout(env,1) # compensate for randomisation time-out

    # pretriage medics
    if paramset.PreTriage
        doctorlastvictim = @process informresourceLastVictim(env, mediclist.onSiteDoctor)
        nurselastvictim = @process informresourceLastVictim(env, mediclist.onSiteNurse)
        @yield doctorlastvictim & nurselastvictim
        newdr = @process addMedic(env, mediclist.ccpDoctor, 1)
        newnu = @process addMedic(env, mediclist.ccpNurse, 1)
        @yield newdr & newnu
        ResourceLogOutputDB(outputDB, "onSiteDoctor", mediclist.onSiteDoctor.capacity, mediclist.onSiteDoctor.level, now(env), "Finished",paramset)
        ResourceLogOutputDB(outputDB, "onSiteNurse", mediclist.onSiteNurse.capacity, mediclist.onSiteNurse.level, now(env), "Finished",paramset)
        ResourceLogOutputDB(outputDB, "ccpDoctor", mediclist.ccpDoctor.capacity, mediclist.ccpDoctor.level, now(env), "Transfer",paramset)
        ResourceLogOutputDB(outputDB, "ccpNurse", mediclist.ccpNurse.capacity, mediclist.ccpNurse.level, now(env), "Transfer",paramset)
    end

    # firefighters

    # @yield timeout(env, SRArrivalUpperBound*(1+OtherTimeBound))

    while firefighters.capacity > 0
        firefighterslastvictim = @process informresourceLastVictim(env, firefighters)
        @yield firefighterslastvictim
        ResourceLogOutputDB(outputDB, "firefighters", firefighters.capacity, firefighters.level, now(env), "Finished",paramset)
    end


    for i=1:consts.nrTriageMMT
        doctorlastvictim = @process informresourceLastVictim(env, mediclist.ccpDoctor)
        nurselastvictim = @process informresourceLastVictim(env, mediclist.ccpNurse)
        @yield doctorlastvictim & nurselastvictim
        newmmt = @process addMedic(env, mediclist.fmpMMT, 1)
        @yield newmmt

        ResourceLogOutputDB(outputDB, "ccpDoctor", mediclist.ccpDoctor.capacity, mediclist.ccpDoctor.level, now(env), "Finished",paramset)
        ResourceLogOutputDB(outputDB, "ccpNurse", mediclist.ccpNurse.capacity, mediclist.ccpNurse.level, now(env), "Finished",paramset)
        ResourceLogOutputDB(outputDB, "fmpMMT", mediclist.fmpMMT.capacity, mediclist.fmpMMT.level, now(env), "Transfer",paramset)

    end


    # with stay and play: ambulances small noria and treatment at fmp
    if policy == "StayPlay" 
        # small noria
        for i in 1:paramset.nrAmbuFW #civilian assets
            ambulastvictim = @process informresourceLastVictim(env, ambulist.ForwardMedevac)
            @yield ambulastvictim
            newamb = @process addMedic(env, ambulist.TacMedevac, 1)
            @yield newamb
            ResourceLogOutputDB(outputDB, "ForwardMedevac", ambulist.ForwardMedevac.capacity, ambulist.ForwardMedevac.level, now(env), "Finished",paramset)
            ResourceLogOutputDB(outputDB, "TacMedevac", ambulist.TacMedevac.capacity, ambulist.TacMedevac.level, now(env), "Transfer",paramset)
        end
        for i in 1:consts.nrAmbuPreliminary #military medical assets
            medevaclastvictim = @process informresourceLastVictim(env, medevac.preliminary)
            @yield medevaclastvictim
            newamb = @process addMedic(env, medevac.primary, 1)
            @yield newamb
            ResourceLogOutputDB(outputDB, "preliminary", medevac.preliminary.capacity, medevac.preliminary.level, now(env), "Finished",paramset)
            ResourceLogOutputDB(outputDB, "primary", medevac.primary.capacity, medevac.primary.level, now(env), "Transfer",paramset)
        end
        # treatment at FMP
        # After last victim no relocation, but all medical personnel may participate in transport
        mediclist.minDoctorFMP = 0
        mediclist.minNurseFMP = 0
    end

     # signal ok to finish transport
    @yield request(mediclist.fmpMMT, priority = 9999)
    @yield release(mediclist.fmpMMT, priority = 9999)
    ambulist.allready = true

end

#DB RELATED FUNCTIONS##

function createOutputDataBase(consts::Constants,baseoutputname::String)
    outputDB = SQLite.DB(baseoutputname)


    # Victim Log
    SQLite.drop!(outputDB, "VictimLog", ifexists = true )
    cmdSQL = "CREATE TABLE VictimLog(
    SimulationRun Int64,
    VictimID Int64,
    Time Float64,
    Action Int64,
    SimScore Float64,
    patientX Float64,
    patientY Float64,
    triage Int64,
    mobility Int64,
    bleed Int64,
    hospid Int64,
    gamma Float64,
    bv Float64)"
    SQLite.execute( outputDB, cmdSQL )

    # Resource Log
    SQLite.drop!(outputDB, "ResourceLog", ifexists = true )
    cmdSQL = "CREATE TABLE ResourceLog(
    SimulationRun int,
    Resource varchar(16),
    Capacity int,
    Level int,
    Time float,
    Action varchar(16) )"
    SQLite.execute( outputDB, cmdSQL )

    # Overview
    SQLite.drop!(outputDB, "Overview", ifexists = true )
    cmdSQL = "CREATE TABLE Overview(
    SimulationRun int,
    Seed int,
    Tourniquet varchar(16),
    Mascal varchar(16),
    Policy varchar(16),
    PreTriage varchar(16),
    Pretriagelocation varchar(16),
    Triagelevel varchar(16),
    Transpsupervision varchar(16),
    Hospdistr varchar(16),
    Hospcapacity varchar(16),
    nrdeathbeforetriage int,
    nrdeathbeforetreattransp int,
    nrdeathbeforehosp int,
    nrdeathinhosp int,
    nrdeathtotal int,
    nrvictimcheck int,
    nrmedicscheck int,
    nrambucheck int
    )  "
    SQLite.execute( outputDB, cmdSQL )

    return outputDB
end
function VictimLogOutputDB(outputDB :: SQLite.DB, VictimID ::Int64,  Time ::Float64, Action ::Int64 , SimScore :: Float64,patientX :: Float64,patientY ::Float64,triage::Int64,mobility::Int64,bleed::Int64,paramset::Parameterset,hospid::Int64,gamma::Float64,bv::Float64)
    cmdSQL =  "INSERT INTO VictimLog values ('$(paramset.simurun)', '$VictimID',  '$Time', '$Action','$SimScore','$patientX','$patientY','$triage','$mobility','$bleed','$hospid','$gamma', '$bv')"
    SQLite.execute( outputDB, cmdSQL )
end
function ResourceLogOutputDB(outputDB :: SQLite.DB, ResourceName ::String, Capacity ::Int64, Level ::Int64, Time ::Float64, Action ::String,paramset::Parameterset)
    cmdSQL =  "INSERT INTO ResourceLog values ('$(paramset.simurun)', '$ResourceName', '$Capacity', '$Level', '$Time', '$Action' )"
    SQLite.execute( outputDB, cmdSQL )
end
function safe_query(db, query)
    df = DataFrame(DBInterface.execute(db, query))
    return df
end
function AddOutputOverviewDB(paramset :: Parameterset, outputDB :: SQLite.DB, mediclist ::MedicList, ambulances ::AmbuList,consts)

    simurun = paramset.simurun
    if consts.logg
        death_query = "SELECT VictimID, Time FROM VictimLog WHERE Action = 188 AND SimulationRun = '$simurun'"
        deathdata = unique(safe_query(outputDB, death_query))
        rename!(deathdata, :Time => :TimeOfDeath)
        ccp_query = "SELECT VictimID, Time FROM VictimLog WHERE Action = 1 AND SimulationRun = '$simurun'"
        ccpdata = safe_query(outputDB, ccp_query)
        if !isempty(deathdata) && !isempty(ccpdata)
            combined_ccp_death = innerjoin(deathdata, ccpdata, on = :VictimID, makeunique=true)
            nrdeathbeforetriage = sum(combined_ccp_death.TimeOfDeath .< combined_ccp_death.Time)
        else
            nrdeathbeforetriage = 0
        end
           
        hosp_query = "SELECT VictimID, Time FROM VictimLog WHERE Action = 17 AND SimulationRun = '$simurun'"
        hospdata = safe_query(outputDB, hosp_query)
        if !isempty(deathdata) && !isempty(hospdata)
            combined_hosp_death = innerjoin(deathdata, hospdata, on = :VictimID, makeunique=true)
            nrdeathbeforehosp = sum(combined_hosp_death.TimeOfDeath .< combined_hosp_death.Time)
            nrdeathinhosp = nrow(combined_hosp_death) - nrdeathbeforehosp
        else
            nrdeathbeforehosp = 0
            nrdeathinhosp = 0
        end
        nrdeathtotal = nrow(deathdata)
        nrdeathbeforetreattransp = nrdeathtotal - nrdeathbeforetriage - nrdeathbeforehosp - nrdeathinhosp
        nrdeathbeforetreattransp = max(nrdeathbeforetreattransp, 0)  # Ensure non-negative value
        home_query = "SELECT VictimID FROM VictimLog WHERE Action = 0 AND SimulationRun = '$simurun'"
        nrhome = nrow(safe_query(outputDB, home_query))
        victim_query = "SELECT VictimID, Time FROM VictimLog WHERE Action = 42 AND SimulationRun = '$simurun'"
        victimdata = safe_query(outputDB, victim_query)
        nrvictimcheck = nrow(victimdata)
    else
        deathdata =  unique(DataFrame(DBInterface.execute( outputDB, "SELECT VictimID, Time FROM VictimLog WHERE Action = 188 AND SimulationRun = '$simurun'" )))
        rename!(deathdata, :Time => :TimeOfDeath)
        nrvictimcheck = 0
        nrdeathbeforetriage = 0
        nrdeathinhosp =0
        nrdeathbeforehosp =0
        nrdeathbeforetreattransp =0
        nrdeathtotal=size(deathdata,1)

    end
        
    nrmedicscheck = mediclist.onSiteDoctor.capacity + mediclist.onSiteNurse.capacity + mediclist.ccpDoctor.capacity + mediclist.ccpNurse.capacity + mediclist.fmpMMT.capacity*2 + mediclist.fmpDoctor.capacity + mediclist.fmpNurse.capacity + mediclist.nrDoctorTransport + mediclist.nrNurseTransport
    nrambucheck = ambulances.ForwardMedevac.capacity + ambulances.TacMedevac.capacity

    ResourceLogOutputDB(outputDB, "onSiteDoctor", mediclist.onSiteDoctor.capacity, mediclist.onSiteDoctor.level, 999.9, "End",paramset)
    ResourceLogOutputDB(outputDB, "onSiteNurse", mediclist.onSiteNurse.capacity, mediclist.onSiteNurse.level, 999.9, "End",paramset)
    ResourceLogOutputDB(outputDB, "ccpDoctor", mediclist.ccpDoctor.capacity, mediclist.ccpDoctor.level, 999.9, "End",paramset)
    ResourceLogOutputDB(outputDB, "ccpNurse", mediclist.ccpNurse.capacity, mediclist.ccpNurse.level, 999.9, "End",paramset)
    ResourceLogOutputDB(outputDB, "fmpMMT", mediclist.fmpMMT.capacity, mediclist.fmpMMT.level, 999.9, "End",paramset)
    ResourceLogOutputDB(outputDB, "fmpDoctor", mediclist.fmpDoctor.capacity, mediclist.fmpDoctor.level, 999.9, "End",paramset)
    ResourceLogOutputDB(outputDB, "fmpNurse", mediclist.fmpNurse.capacity, mediclist.fmpNurse.level, 999.9, "End",paramset)
    ResourceLogOutputDB(outputDB, "nrDoctorTransport", mediclist.nrDoctorTransport, mediclist.nrDoctorTransport, 999.9, "End",paramset)
    ResourceLogOutputDB(outputDB, "nrNurseTransport", mediclist.nrNurseTransport, mediclist.nrNurseTransport, 999.9, "End",paramset)

    data = [simurun, paramset.seed, paramset.Tourniquet,paramset.mascal, paramset.Policy, paramset.PreTriage, paramset.PreTriageLocation, paramset.TriageLevel, paramset.TranspsupervisionLevel, paramset.HospitalDistribution, paramset.HospitalCapacity, nrdeathbeforetriage, nrdeathbeforetreattransp, nrdeathbeforehosp, nrdeathinhosp, nrdeathtotal, nrvictimcheck, nrmedicscheck, nrambucheck]


    datastring = "'$(data[1])'"
        for i in data[2:end]
            datastring *= ", '$i'"
        end
        cmdSQL =  "INSERT INTO Overview values ("*datastring*")"
        SQLite.execute( outputDB, cmdSQL )
end

#SCHEDULING FUNCTIONS

@resumable function SRArrival(env :: Environment, firefighters :: Resource, paramset :: Parameterset,arrivalTime::Float64,consts::Constants,outputDB)
   
    @yield timeout(env, timerandomised(arrivalTime, "Other",consts))
    
    ResourceLogOutputDB(outputDB, "Firefighter_team", firefighters.capacity, firefighters.level, now(env), "Arrival",paramset)
   
    
    
end
@resumable function scheduleSR(env :: Environment, firefighters :: Resource, paramset :: Parameterset,consts::Constants,outputDB,inputDB,disasters,m,threat,alertTimes,offroad::Bool)

    
 
    dfAMBselection = getTimeline(inputDB,disasters,m,"Fire",consts,paramset,alertTimes,offroad)
    threats=[]
    for i in 1:length(disasters)
        push!(threats, "threat$i")
    end

    # Iterate over threats
    for (threat_index, threat_name) in enumerate(threats)
        if threat == threat_index - 1
            for i in eachindex(dfAMBselection)
                if dfAMBselection[i][6] == threat_name # Firefighters to the corresponding threat site
                    alertTime = timerandomised(1.0, "Other", consts)
                    travelTime = dfAMBselection[i][1]
                    arrivalTime = alertTime + timerandomised(travelTime, "Travel", consts)
                    @process SRArrival(env, firefighters, paramset, arrivalTime, consts, outputDB)
                end
            end
        end
    end
   
end
@resumable function scheduleMMT(env :: Environment, mediclist :: MedicList, ambulist :: AmbuList, paramset :: Parameterset, inputDB :: SQLite.DB,consts::Constants,outputDB)

    MMTdata = DataFrame(CSV.File("timeline.txt"; normalizenames=true))
   
    for i in 1:size(MMTdata,1)
        if MMTdata[i,5] == 1
        
            alertTime = timerandomised(1.0, "Other",consts)
            
            mmtName = MMTdata[i,4]
            
            
            travelTime = MMTdata[i,3]
            arrivalTime = alertTime + timerandomised(travelTime, "Travel",consts)
            role =""
            @process mmtArrival(env, mediclist, ambulist, paramset, arrivalTime,role,consts,outputDB)
        end

    end
end
@resumable function scheduleMEDEVAC(env :: Environment, mediclist::MedicList, medevac :: Medevac, paramset :: Parameterset, inputDB :: SQLite.DB,consts::Constants,outputDB)

    dfAMBselection = DataFrame(CSV.File("timelineMil.txt"; normalizenames=true))
    
    for i in 1:size(dfAMBselection,1)
        if dfAMBselection[i,5] == 2 #armored
            ambName = dfAMBselection[i,4]
        
            alertTime = timerandomised(1.0, "Other",consts)
            travelTime = dfAMBselection[i,3]
            arrivalTime = alertTime + timerandomised(travelTime, "Travel",consts)
                  
            @process milArrival(env, mediclist, medevac, paramset, arrivalTime, "armored",consts,outputDB)
        end
        # nr units
      
    end

end
function getTimeline(inputDB::SQLite.DB,disasters::Array,m,unit::String,consts,paramset,alertTimes,offroad::Bool)
    #read MTF coords
    timeline=[]
    timeline1=[]
    timeline2=[]
    sorted_timeline =[]
    timelines = []
    threats =[]
  
    if !offroad
  
        if unit == "Amb"
            MTFlist = DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Ambulance" ))
            thresholdTime = 8000.0
            if length(disasters) == 1

                #push!(timeline,["driveTime"   "distance"   "Name" "id" "type"])
                for i in 1:min(paramset.nrAmbuFW,size(MTFlist,1))
                    originECF = ECEFfromLLA(wgs84)(Geodesy.LLA(disasters[1][1],disasters[1][2]))
                    hospECF = ECEFfromLLA(wgs84)(Geodesy.LLA(MTFlist[i,:].LAT,MTFlist[i,:].LON))
                
                    distToHosp=OpenStreetMapX.distance(originECF.x,originECF.y,originECF.z,hospECF.x,hospECF.y,hospECF.z)
              
                    # if paramset.simurun ==1
                    #     push!(routes,route)
                    # end OFFROAD -> no routes
                    driveTime = (distToHosp*1.5/(MTFlist[i,:].maxspeed*(16.6666))) #time in mins is distance in m divided by speed in m/min
                   
                
                    if driveTime < thresholdTime
                        push!(timeline,[driveTime   distToHosp   MTFlist[i,:].name i 2    "threat1"])
                    end
                end
                # Get the indices that would sort the first column
                sorted_indices = sortperm([row[1] for row in timeline])

                # Rearrange the timeline based on the sorted indices
                sorted_timeline = timeline[sorted_indices]
                # Get the indices that would sort the first column
                
                writedlm("timelineEMS.txt",sorted_timeline)
            # push!(timeline1,["driveTime"   "distance"   "Name" "id" "type"])
                #push!(timeline,["driveTime"   "distance"   "Name" "id" "type"])
                # for i in 1:min(paramset.nrAmbuFW,size(MTFlist,1))
            
                #     hospNode = point_to_nodes((MTFlist[i,:].LAT,MTFlist[i,:].LON),m)
                #     originNode = point_to_nodes((disasters[1][1],disasters[1][2]),m)
                #     route, distance,Time = OpenStreetMapX.fastest_route(m,originNode,hospNode)  
                #     driveTime = (Time/60.0)*1.5
                
                #     push!(timeline,[driveTime   distance   MTFlist[i,:].Name i 2    "threat1"])
                # end
                # # Get the indices that would sort the first column
                # sorted_indices = sortperm([row[1] for row in timeline])

                # # Rearrange the timeline based on the sorted indices
                # sorted_timeline = timeline[sorted_indices]
                # Get the indices that would sort the first column
                return sorted_timeline
            end
        # push!(timeline1,["driveTime"   "distance"   "Name" "id" "type"])
            if length(disasters) >= 2
                MTFlist = DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Hospital" ))
                for k in 1:length(disasters)
                    for i in 1:min(paramset.nrAmbuFW,size(MTFlist,1))
                
                        hospNode = point_to_nodes((MTFlist[i,:].LAT,MTFlist[i,:].LON),m)
                        originNode = point_to_nodes((disasters[k][1],disasters[k][2]),m)
                        route, distance,Time = OpenStreetMapX.fastest_route(m,originNode,hospNode)  
                        if Time ==0                                                     
                            hospECF = ECEFfromLLA(wgs84)(Geodesy.LLA(MTFlist[i, :].LAT, MTFlist[i, :].LON))
                            originECF = ECEFfromLLA(wgs84)(Geodesy.LLA(disasters[k][1], disasters[k][2]))
                            distanceA = OpenStreetMapX.distance(OpenStreetMapX.ECEF(hospECF.x, hospECF.y,hospECF.z), OpenStreetMapX.ECEF(originECF.x, originECF.y,originECF.z))
                            driveTime = (distanceA * 2.0)/(30.0*(16.6666)) # time in minutes
                            println("$distanceA, $driveTime $hospECF $originECF")
                                
                          

                        else
                            driveTime = (Time/60.0)*1.5 + alertTimes[k]
                        end
                    
                        push!(timeline1,[driveTime distance MTFlist[i,:].Name   i   k  "threat$k"])
                    end
                    sorted_indices1 = sortperm([row[1] for row in timeline1])
                    sorted_timeline1 = timeline1[sorted_indices1]

                    # for i in 1:min(paramset.nrAmbuFW,size(MTFlist,1))
                    
                    #     hospNode = point_to_nodes((MTFlist[i,:].LAT,MTFlist[i,:].LON),m)
                    #     originNode = point_to_nodes((disasters[3]),m)
                    #     route, distance,Time = OpenStreetMapX.fastest_route(m,originNode,hospNode)  
                    #     driveTime = (Time/60.0)*1.5
                    
                    #     push!(timeline2,[driveTime distance MTFlist[i,:].Name   i   2  "threat3"])
                    # end
                    # sorted_indices2 = sortperm([row[1] for row in timeline2])
                    # sorted_timeline2 = timeline2[sorted_indices2]


                    #writedlm("timeline2.txt",sorted_timeline1)
                    timelineall = []
                    timelineall = append!(timelineall,sorted_timeline,sorted_timeline1)
                    #timelineall = append!(timelineall,sorted_timeline,sorted_timeline1,sorted_timeline2)
                    writedlm("timelineEMS.txt",timelineall)
                    end
                return timelineall
            end
            if length(disasters) == 1
                return sorted_timeline

            end
        end
        if unit == "Fire"
            FFlist = DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Firefighters" ))
            
            if length(disasters) == 1
                for i in 1:size(FFlist,1)
            
                    hospNode = point_to_nodes((FFlist[i,:].LAT,FFlist[i,:].LON),m)
                    originNode = point_to_nodes((disasters[1][1],disasters[1][2]),m)
                    route, distance,Time = OpenStreetMapX.fastest_route(m,originNode,hospNode)  
                    driveTime = (Time/60.0)*1.5
                
                    push!(timeline,[driveTime   distance   FFlist[i,:].Name i 1    "threat1"])
                end
                # Get the indices that would sort the first column
                sorted_indices = sortperm([row[1] for row in timeline])

                # Rearrange the timeline based on the sorted indices
                sorted_timeline = timeline[sorted_indices]
                # Get the indices that would sort the first column
                
                writedlm("timeline1ff.txt",sorted_timeline)
      
            end
            if length(disasters) > 2
                # Define threats explicitly
                for i in 1:10
                    push!(threats, "threat$i")
                end

                # Loop through each threat
                for (threat_index, threat_name) in enumerate(threats)
                    if threat_index + 1 <= length(disasters) # Ensure there are enough disasters defined
                        timeline = []
                        originNode = point_to_nodes((disasters[threat_index + 1][1],disasters[threat_index + 1][2]), m) # +1 because threats start from threat2

                        # Loop through each firefighter in FFlist
                        for i in 1:size(FFlist, 1)
                            hospNode = point_to_nodes((FFlist[i, :].LAT, FFlist[i, :].LON), m)
                            route, distance, Time = OpenStreetMapX.fastest_route(m, originNode, hospNode)
                            driveTime = (Time / 60.0) * 1.5
                            push!(timeline, [driveTime distance FFlist[i, :].Name i 1 threat_name])
                        end

                        # Sort the timeline based on driveTime
                        sorted_indices = sortperm([row[1] for row in timeline])
                        sorted_timeline = timeline[sorted_indices]
                        push!(timelines, sorted_timeline)
                    end
                end

                # Combine all sorted timelines
                timelineall = vcat(timelines...)

                # Write the combined timeline to a file
                writedlm("timelineFF.txt", timelineall)
            #     for i in 1:size(FFlist,1)
            #         # println("hosp $i : $(MTFlist[i].LAT)  $(MTFlist[i].LON)")
            #         hospNode = point_to_nodes((FFlist[i,:].LAT,FFlist[i,:].LON),m)
            #         originNode = point_to_nodes((disasters[2]),m)
            #         route, distance,Time = OpenStreetMapX.fastest_route(m,originNode,hospNode)  
            #         driveTime = (Time/60.0)*1.5
            #         #println("MTF$i $driveTime to zaventem")
            #         push!(timeline1,[driveTime distance FFlist[i,:].Name   i   1  "threat2"])
            #     end
            #     sorted_indices1 = sortperm([row[1] for row in timeline1])
            #     sorted_timeline1 = timeline1[sorted_indices1]

            #     for i in 1:size(FFlist,1)
            #         # println("hosp $i : $(MTFlist[i].LAT)  $(MTFlist[i].LON)")
            #         hospNode = point_to_nodes((FFlist[i,:].LAT,FFlist[i,:].LON),m)
            #         originNode = point_to_nodes((disasters[3]),m)
            #         route, distance,Time = OpenStreetMapX.fastest_route(m,originNode,hospNode)  
            #         driveTime = (Time/60.0)*1.5
            #         #println("MTF$i $driveTime to zaventem")
            #         push!(timeline1,[driveTime distance FFlist[i,:].Name   i   1  "threat3"])
            #     end
            #     sorted_indices2 = sortperm([row[1] for row in timeline2])
            #     sorted_timeline2 = timeline2[sorted_indices2]
        

            # #  writedlm("timeline2ff.txt",sorted_timeline1)
            #     timelineall = []
            #     timelineall = append!(timelineall,sorted_timeline,sorted_timeline1,sorted_timeline2)
            # # writedlm("timelineFF.txt",timelineall)
                 return timelineall
            end
            if length(disasters) == 1
                return sorted_timeline
            end

        end
    else
        if unit == "Amb"
            MTFlist = DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Ambulance" ))
            thresholdTime = 8000.0
            if length(disasters) == 1

                #push!(timeline,["driveTime"   "distance"   "Name" "id" "type"])
                for i in 1:min(paramset.nrAmbuFW,size(MTFlist,1))
                    originECF = ECEFfromLLA(wgs84)(Geodesy.LLA(disasters[1][1],disasters[1][2]))
                    hospECF = ECEFfromLLA(wgs84)(Geodesy.LLA(MTFlist[i,:].LAT,MTFlist[i,:].LON))
                
                    distToHosp=OpenStreetMapX.distance(originECF.x,originECF.y,originECF.z,hospECF.x,hospECF.y,hospECF.z)
              
                    # if paramset.simurun ==1
                    #     push!(routes,route)
                    # end OFFROAD -> no routes
                    driveTime = (distToHosp*2/(50*(16.6666))) #time in mins is distance in m divided by speed in m/min times a correction factor depending on traffic (to calibrate)
                    
                
                    if driveTime < thresholdTime
                        push!(timeline,[driveTime   distToHosp   MTFlist[i,:].Name i 2    "threat1"])
                    end
                end
                # Get the indices that would sort the first column
                sorted_indices = sortperm([row[1] for row in timeline])

                # Rearrange the timeline based on the sorted indices
                sorted_timeline = timeline[sorted_indices]
                # Get the indices that would sort the first column
                
                writedlm("timelineEMS.txt",sorted_timeline)
            # push!(timeline1,["driveTime"   "distance"   "Name" "id" "type"])
            end
            
            if length(disasters) >= 2
                # Define threats explicitly
                
               

                for i in 1:length(disasters)
                    push!(threats, "threat$i")
                end

                # Loop through each threat
                for (threat_index, threat_name) in enumerate(threats)
                    if threat_index <= length(disasters)
                        timeline = []
                        originECF = ECEFfromLLA(wgs84)(Geodesy.LLA(disasters[threat_index][1], disasters[threat_index][2]))

                        # Loop through each ambulance in MTF list
                        for i in 1:min(paramset.nrAmbuFW, size(MTFlist, 1))
                            hospECF = ECEFfromLLA(wgs84)(Geodesy.LLA(MTFlist[i, :].LAT, MTFlist[i, :].LON))

                            distToHosp = OpenStreetMapX.distance(originECF.x, originECF.y, originECF.z, hospECF.x, hospECF.y, hospECF.z)
                            driveTime = (distToHosp * 2.0)/(MTFlist[i, :].maxspeed*(16.6666)) # time in minutes

                            if driveTime < thresholdTime
                                push!(timeline, [driveTime distToHosp MTFlist[i, :].name i 2 threat_name originECF.x originECF.y])
                            end
                        end

                        # Sort the timeline based on driveTime
                        sorted_indices = sortperm([row[1] for row in timeline])
                        sorted_timeline = timeline[sorted_indices]

                        push!(timelines, sorted_timeline)
                    end
                end

                # Combine all sorted timelines
               
                timelineall = vcat(timelines...)

                # Write the combined timeline to a file
                writedlm("timelineEMS.txt", timelineall)

                # Return the combined timeline
                return timelineall
                
            end
            if length(disasters) == 1
                return sorted_timeline

            end
        end
        if unit == "Fire"
            FFlist = DataFrame(DBInterface.execute( inputDB, "SELECT * FROM Firefighters" ))
            timelines=[]
            if length(disasters) == 1
                for i in 1:size(FFlist,1)
                    
                    originECF = ECEFfromLLA(wgs84)(Geodesy.LLA(disasters[1][1],disasters[1][2]))
                    hospECF = ECEFfromLLA(wgs84)(Geodesy.LLA(FFlist[i,:].LAT,FFlist[i,:].LON))
                
                    distToHosp=OpenStreetMapX.distance(originECF.x,originECF.y,originECF.z,hospECF.x,hospECF.y,hospECF.z)
                
                    driveTime = (distToHosp*2.0/(30*(16.6666))) #time in mins is distance in m divided by speed in m/min
            
                    
              
                    push!(timeline,[driveTime   distToHosp   FFlist[i,:].Name i 1    "threat1"])
                end
                # Get the indices that would sort the first column
                sorted_indices = sortperm([row[1] for row in timeline])

                # Rearrange the timeline based on the sorted indices
                sorted_timeline = timeline[sorted_indices]
                # Get the indices that would sort the first column
                
                writedlm("timeline1ff.txt",sorted_timeline)
         
            end
            if length(disasters) > 2
                for i in 1:10
                    push!(threats, "threat$i")
                end

                # Constants
                speed_m_per_min = 50 * 16.6666 # speed in m/min 

                # Loop through each threat
                for (threat_index, threat_name) in enumerate(threats)
                    if threat_index + 1 <= length(disasters) # Ensure there are enough disasters defined
                        timeline = []
                        originECF = ECEFfromLLA(wgs84)(Geodesy.LLA(disasters[threat_index + 1][1], disasters[threat_index + 1][2]))

                        # Loop through each firefighter in FFlist
                        #for i in 1:size(FFlist, 1)
                            hospECF = ECEFfromLLA(wgs84)(Geodesy.LLA(FFlist[1, :].LAT, FFlist[1, :].LON))
                            distToHosp = OpenStreetMapX.distance(originECF.x, originECF.y, originECF.z, hospECF.x, hospECF.y, hospECF.z)
                            driveTime = (distToHosp * 2.0) / speed_m_per_min # time in minutes
                            push!(timeline, [driveTime distToHosp FFlist[1, :].Name 1 1 threat_name])
                        #end

                        # Sort the timeline based on driveTime
                        sorted_indices = sortperm([row[1] for row in timeline])
                        sorted_timeline = timeline[sorted_indices]
                        push!(timelines, sorted_timeline)
                    end
                end

                # Combine all sorted timelines
                timelineall = vcat(timelines...)

                # Write the combined timeline to a file
                writedlm("timelineFF.txt", timelineall)

                # Return the combined timeline
                return timelineall
           
            end
            if length(disasters) == 1
                return sorted_timeline
            end

        end

    end
        
end
@resumable function scheduleBOAT(env :: Environment, boatlist :: BoatList, paramset :: Parameterset, inputDB :: SQLite.DB,consts::Constants,outputDB,disasters,m,threat)

    alertTime = timerandomised(2.0, "Other",consts)
    arrivalTime = alertTime 
    @process boatArrival(env, boatlist, paramset, arrivalTime, consts,outputDB)

end
@resumable function scheduleAMB(env :: Environment, ambulist :: AmbuList,mediclist::MedicList, paramset :: Parameterset, inputDB :: SQLite.DB,consts::Constants,outputDB,disasters,m,threat,alertTimes::Array,offroad::Bool)
    
    dfAMBselection=getTimeline(inputDB,disasters,m,"Amb",consts,paramset,alertTimes,offroad)
   
    
    threats=[]
    for i in 1:length(disasters)
        push!(threats, "threat$i")
    end
    
 # Iterate over threats
    for k in 1:length(threats)
       
        for i in eachindex(dfAMBselection)
            
            if dfAMBselection[i][6] == threats[k] && i < paramset.nrAmbuFW # Ambulance to the corresponding threat site and enough ambulances
                alertTime = timerandomised(alertTimes[k], "Other", consts)
                travelTime = dfAMBselection[i][1]
                arrivalTime = alertTime + timerandomised(travelTime, "Travel", consts)
                
                @process ambuArrival(env, ambulist,mediclist, paramset, arrivalTime, "BLS", consts, outputDB)
            end
        end
       
    end
    
end

#TRANSPORT RESOURCES FUNCTIONS
@resumable function mmtArrival(env ::Environment, mediclist ::MedicList, ambulist :: AmbuList, paramset ::Parameterset, delaytime ::Float64, role ::String,consts::Constants,outputDB)

    @yield timeout(env, delaytime)

    # create MMT
    arrivaltime = now(env)

    # DirMed
    plannedDirMed = 0

    if (plannedDirMed > 0 && role == "DirMed") || ( role == "" && plannedDirMed == 0 && mediclist.dirMed.capacity == 0)
        mediclist.dirMed.capacity += 1
        ResourceLogOutputDB(outputDB, "DirMed", mediclist.dirMed.capacity, mediclist.dirMed.level, now(env), "Arrival",paramset)

    # First Triage, then NUCA
    elseif mediclist.nrTriageMMTArrived < consts.nrTriageMMT
        mediclist.nrTriageMMTArrived += 1

        if paramset.PreTriage && (mediclist.nrTriageMMTArrived == (consts.nrTriageMMT - 1)) # let first triage mmts go to ccp, let last go to site fro pretriage
            # CCP triage start before pre-triage

            mediclist.onSiteDoctor.capacity += 1
            @yield request(mediclist.onSiteDoctor, priority = 1)
            @yield release(mediclist.onSiteDoctor)
            ResourceLogOutputDB(outputDB, "onSiteDoctor", mediclist.onSiteDoctor.capacity, mediclist.onSiteDoctor.level, now(env), "Arrival",paramset)

            mediclist.onSiteNurse.capacity += 1
            @yield request(mediclist.onSiteNurse, priority = 1)
            @yield release(mediclist.onSiteNurse)
            ResourceLogOutputDB(outputDB, "onSiteNurse", mediclist.onSiteNurse.capacity, mediclist.onSiteNurse.level, now(env), "Arrival",paramset)
           # println("On site doctor +1 nurse arrived at $(now(env))!")

        else
           # println("ccp doctor +1 nurse arrived at $(now(env))!")
            mediclist.ccpDoctor.capacity += 1
            @yield request(mediclist.ccpDoctor, priority = 1)
            @yield release(mediclist.ccpDoctor)
            ResourceLogOutputDB(outputDB, "ccpDoctor", mediclist.ccpDoctor.capacity, mediclist.ccpDoctor.level, now(env), "Arrival",paramset)
            mediclist.ccpNurse.capacity += 1
            @yield request(mediclist.ccpNurse, priority = 1)
            @yield release(mediclist.ccpNurse)
            ResourceLogOutputDB(outputDB, "ccpNurse", mediclist.ccpNurse.capacity, mediclist.ccpNurse.level, now(env), "Arrival",paramset)

        end

    # FMP
    else
        if paramset.Policy == "StayPlay" && paramset.buildFMP#FMP exists already
            # Building FMP
            if mediclist.BuildingFMP == false
                mediclist.BuildingFMP = true
                @process startBuildingFMP(env, mediclist,consts)
            end
            # Drive to FMP
            @yield timeout(env, timerandomised(consts.DriveTimeFMP, "Other",consts))
        end
       # println("fmp doctor +1 nurse arrived at $(now(env))!")
        mediclist.fmpMMT.capacity += 1
        @yield request(mediclist.fmpMMT, priority = 1)
        @yield release(mediclist.fmpMMT)
        ResourceLogOutputDB(outputDB, "fmpMMT", mediclist.fmpMMT.capacity, mediclist.fmpMMT.level, now(env), "Arrival",paramset)

        @yield release(ambulist.transportchecktrigger)

    end
end
@resumable function milArrival(env :: Environment, mediclist::MedicList,medevac :: Medevac, paramset :: Parameterset, delaytime :: Float64, ambutype :: String,consts :: Constants,outputDB)

    @yield timeout(env, delaytime)

    if ambutype == "armored"
        # create armored ambulance + EMT
        
        if (paramset.Policy == "StayPlay" ) && medevac.preliminary.capacity < consts.nrAmbuPreliminary

            medevac.preliminary.capacity += 1
            @yield request(medevac.preliminary, priority = 1)
            @yield release(medevac.preliminary)
            ResourceLogOutputDB(outputDB, "preliminary",medevac.preliminary.capacity, medevac.preliminary.level, now(env), "Arrival",paramset)
        else
         
            medevac.primary.capacity += 1
            @yield request(medevac.primary, priority = 1)
            @yield release(medevac.primary)
            ResourceLogOutputDB(outputDB, "primary",medevac.primary.capacity, medevac.primary.level, now(env), "Arrival",paramset)
            @yield release(medevac.transportchecktrigger)
        end
    end
end
@resumable function boatArrival(env :: Environment, boatList :: BoatList, paramset :: Parameterset, delaytime :: Float64, consts::Constants,outputDB)

    @yield timeout(env, delaytime)

            # create zodiac
        
    for i in 1:consts.nrBoats    
              
         
            
     
            boatList.ForwardMedevac.capacity += 1
            # @yield request(boatList.ForwardMedevac, priority = 1)
            # @yield release(boatList.ForwardMedevac)
            ResourceLogOutputDB(outputDB, "zodiac",boatList.ForwardMedevac.capacity, boatList.ForwardMedevac.level, now(env), "Arrival",paramset)
   
    end
    
    
end
@resumable function ambuArrival(env :: Environment, ambulist :: AmbuList,mediclist::MedicList, paramset :: Parameterset, delaytime :: Float64, ambutype :: String,consts::Constants,outputDB)

    @yield timeout(env, delaytime)

    if ambutype == "BLS"
        # create BLS ambulance + 2 EMT
        
        if (paramset.Policy == "StayPlay" ) 
    
           # @yield request(ambulist.ForwardMedevac, priority = 1)
          #  @yield release(ambulist.ForwardMedevac)
            ResourceLogOutputDB(outputDB, "ForwardMedevac",ambulist.ForwardMedevac.capacity, ambulist.ForwardMedevac.level, now(env), "Arrival",paramset)
            @process addMedic(env,mediclist.ccpNurse,1)
        else
          
            
            # @yield request(ambulist.TacMedevac, priority = 1)
            # @yield release(ambulist.TacMedevac)
            ResourceLogOutputDB(outputDB, "TacMedevac",ambulist.TacMedevac.capacity, ambulist.TacMedevac.level, now(env), "Arrival",paramset)
            @process addMedic(env,mediclist.ccpNurse,1)
            # @yield release(ambulist.transportchecktrigger)
        end
    end
end

end