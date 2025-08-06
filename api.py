from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel
from fastapi import Body
from typing import List, Optional, Dict, Any
import subprocess
import json
from datetime import datetime
from pathlib import Path
import glob
import os
import sqlite3
from fastapi import Query
from fastapi import Form
import logging


app = FastAPI()

# Middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)
class ThreatConfig(BaseModel):
    id: str
    hits: int
    victimsPerHit: int
    radius: float
    duration: float
    scheduledTime: float

class ScenarioRequest(BaseModel):
    scenario: str
    threatConfig: Optional[ThreatConfig] = None
# Static and template files
BASE_DIR = Path("/home/pims/SimedisAPI") #change to local API directory
app.mount("/icons", StaticFiles(directory="icons"), name="icons")
templates = Jinja2Templates(directory=BASE_DIR / "templates")

# Serve HTML
@app.get("/", response_class=HTMLResponse)
async def get_index(request: Request):
    return templates.TemplateResponse("index.html", {
        "request": request,
        "now": datetime.utcnow()  # pass current datetime to template
    })

@app.get("/db_browser", response_class=HTMLResponse)
async def db_browser(request: Request):
    now = datetime.now()
    folder_path = "/home/pims/SimedisAPI/"
    db_files = [
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.endswith(".sqlite")
    ]
    return templates.TemplateResponse("db_browser.html", {
        "request": request,
        "db_files": db_files,
        "folder_path": folder_path,   
        "selected_db": None,
        "tables": [],
        "rows": [],
        "headers": [],
        "selected_table": None,
        "now": datetime.now()
    })
# API to run Julia script
@app.post("/api/threatconfig")
async def apply_threat_config(request: Request):
    data = await request.json()
    threat_id = str(data.get("id"))

    if not threat_id:
        return JSONResponse(status_code=400, content={"error": "Missing threat ID"})

    # Load existing config if exists
    config_path = "threatconfig.json"
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            all_configs = json.load(f)
    else:
        all_configs = {}

    # Update/insert the threat config
    all_configs[threat_id] = {
        "scheduledTime": data["scheduledTime"],
        "numberOfHits": data["hits"],
        "victimsPerHit": data["victimsPerHit"],
        "radius": data["radius"]
    }

    # Save all threat configs back
    with open(config_path, "w") as f:
        json.dump(all_configs, f, indent=2)

    return JSONResponse(content={"status": "success", "message": f"Saved config for threat {threat_id}"})
@app.post("/run_scenario")
async def run_scenario(data: dict = Body(...)):
    scenario_name = data.get("scenario", "demo")
    threat_config = data.get("threatConfig")

    julia_script = f"/home/pims/SimedisAPI/{scenario_name}.jl"
    if not os.path.exists(julia_script):
        return JSONResponse(status_code=404, content={"error": f"{julia_script} not found."})

    
    if threat_config:
        with open("threatConfig.json", "w") as f:
            json.dump(threat_config, f, indent=2)
        print(f"[INFO] Saved threat config: {threat_config}")
    else:
        print("[INFO] No threatConfig provided")

    try:
        result = subprocess.run(
            ["julia", julia_script],
            capture_output=True,
            text=True
        )
        return {
            "stdout": result.stdout,
            "stderr": result.stderr,
            "success": result.returncode == 0
        }

    except Exception as e:
        return {
            "stdout": "",
            "stderr": str(e),
            "success": False
        }

@app.post("/generate_patients")
async def generate_patients_api():
    try:
        result = subprocess.run(
            ["julia", "/home/pims/SimedisAPI/genPatients.jl"],
            capture_output=True,
            text=True
        )
        return {
            "stdout": result.stdout,
            "stderr": result.stderr,
            "success": result.returncode == 0
        }
    except Exception as e:
        return {
            "stdout": "",
            "stderr": str(e),
            "success": False
        }
@app.get("/get_victims")
async def get_victims():
    victims = []
    victim_files = glob.glob("/home/pims/SimedisAPI/victims*.txt")
    for file_path in victim_files:
        try:
            with open(file_path, "r") as f:
                lines = f.readlines()
            if not lines:
                continue

            headers = {name: idx for idx, name in enumerate(lines[0].strip().split("\t"))}

            # Defensive check for expected columns
            if "viclat" not in headers or "viclon" not in headers:
                print(f"[WARN] File {file_path} missing 'viclat' or 'viclon'. Headers found: {headers}")
                continue

            for line in lines[1:]:
                parts = line.strip().split("\t")
                if len(parts) < len(headers):
                    continue
                try:
                    lat = float(parts[headers["viclat"]])
                    lng = float(parts[headers["viclon"]])
                    vic_id = parts[headers["victag"]]
                    # Defensive check for valid coordinates
                    if not (-90 <= lat <= 90 and -180 <= lng <= 180):
                        print(f"[WARN] Skipping victim with out-of-range coordinates: {lat}, {lng}")
                        continue

                    triage = parts[headers["triage"]] if "triage" in headers else ""
                    iss = parts[headers["ISS"]] if "ISS" in headers else ""
                    victims.append({
                        "lat": lat,
                        "lng": lng,
                        "id": vic_id,
                        "triage": triage,
                        "ISS": iss
                    })
                except Exception as e:
                    print(f"[ERROR] Parsing line in {file_path}: {line}\nException: {e}")
        except Exception as e:
            print(f"[ERROR] Reading file {file_path}: {e}")

    return {"victims": victims}


def find_sqlite_files(directory: str):
    return [f for f in os.listdir(directory) if f.endswith(".sqlite")]
# Scenario location model
class Location(BaseModel):
    category: str  # threat, hospital, ccp
    subtype: Optional[str] = None
    lat: float
    lng: float

class LocationData(BaseModel):
    locations: List[Location]
    route: Optional[Dict[str, Any]] = None

# Save JSON scenario file

def get_coords(loc):
    if "coords" in loc:
        return loc["coords"]
    elif "coordinates" in loc:
        return loc["coordinates"]
    elif "lat" in loc and "lng" in loc:
        return [loc["lat"], loc["lng"]]
    else:
        print("Bad entry without coordinates:", loc)
        raise ValueError("Missing coordinates in entry")

@app.post("/save_coordinates")
async def save_coordinates(request: Request):
    data = await request.json()
   
    saved_data = {
        "threats": [],
        "ccps": [],
        "hospitals": [],
        "home": None
    }

    home_found = False

    for loc in data.get("locations", []):  # Make sure to access the correct key
        try:
            coords = get_coords(loc)
        except ValueError as e:
            return {"error": str(e)}

        category = loc.get("category", "").lower()
        subtype = loc.get("subtype", "")

        if category == "threat":
            saved_data["threats"].append({"coords": coords, "threatType": subtype})
        elif category == "ccp":
            saved_data["ccps"].append({"coords": coords})
        elif category == "hospital":
            saved_data["hospitals"].append({"coords": coords, "hospitalType": subtype})
        elif category == "home":
            saved_data["home"] = coords
            home_found = True

    if not home_found:
        saved_data["home"] = [50.8467, 4.3525]  # default home coords

    file_path = "scenario.json"
    try:
        with open(file_path, "w") as f:
            json.dump(saved_data, f, indent=2)
        print(f"File saved successfully at {file_path}")
    except Exception as e:
        print(f"Failed to save file: {e}")
        return {"error": "Failed to save file"}

    return {"file": file_path}
@app.post("/db_browser", response_class=HTMLResponse)
async def db_browser_post(
    request: Request,
    db_path: str = Form(...),
    table_name: str = Form(None)
):
    folder_path = "/home/pims/simedis-ui/runSimedis"
    db_files = [
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.endswith(".sqlite")
    ]

    tables = []
    rows = []
    headers = []

    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [row[0] for row in cursor.fetchall()]

        if table_name:
            cursor.execute(f"SELECT * FROM {table_name} LIMIT 100")
            rows = cursor.fetchall()
            headers = [description[0] for description in cursor.description]
        conn.close()
    except Exception as e:
        return HTMLResponse(f"<h3>Error: {e}</h3>")

    return templates.TemplateResponse("db_browser.html", {
        "request": request,
        "db_files": db_files,
        "folder_path": folder_path,
        "selected_db": db_path,
        "tables": tables,
        "rows": rows,
        "headers": headers,
        "selected_table": table_name,
        "now": datetime.now()
    })


@app.get("/sqlite_tables")
async def sqlite_tables(file_path: str = Query(...)):
    target_file = (BASE_DIR / file_path).resolve()
    if not str(target_file).startswith(str(BASE_DIR)) or not target_file.is_file():
        return JSONResponse(status_code=400, content={"error": "Invalid file path"})

    try:
        conn = sqlite3.connect(str(target_file))
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [row[0] for row in cursor.fetchall()]
        conn.close()
        return {"tables": tables}
    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})

@app.get("/sqlite_table_data")
async def sqlite_table_data(file_path: str = Query(...), table_name: str = Query(...), limit: int = 100):
    target_file = (BASE_DIR / file_path).resolve()
    if not str(target_file).startswith(str(BASE_DIR)) or not target_file.is_file():
        return JSONResponse(status_code=400, content={"error": "Invalid file path"})

    try:
        conn = sqlite3.connect(str(target_file))
        cursor = conn.cursor()
        cursor.execute(f"SELECT * FROM {table_name} LIMIT {limit}")
        rows = cursor.fetchall()
        columns = [description[0] for description in cursor.description]
        conn.close()
        return {"columns": columns, "rows": rows}
    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})

@app.get("/browse_files")
async def browse_files(dir_path: str = ""):
    # Sanitize path so no path traversal outside BASE_DIR
    target_path = (BASE_DIR / dir_path).resolve()
    if not str(target_path).startswith(str(BASE_DIR)):
        return JSONResponse(status_code=400, content={"error": "Invalid directory"})

    if not target_path.exists() or not target_path.is_dir():
        return JSONResponse(status_code=404, content={"error": "Directory not found"})

    entries = []
    for entry in target_path.iterdir():
        entries.append({
            "name": entry.name,
            "is_dir": entry.is_dir()
        })

    return {"current_path": str(target_path.relative_to(BASE_DIR)), "entries": entries}