import os
import requests
import earthaccess
import pandas as pd
import datetime as dt
import geopandas
from shapely.geometry import MultiPolygon, Polygon, box
import keyring

# Params =>
lat, lon = 34.0, -118.25
# Temporal Bound - Year, month, day. Hour, minutes, and seconds (ZULU) can also be included 
start_date = dt.datetime(2022, 9, 3)
end_date = dt.datetime(2024, 3, 9, 23, 23, 59)

# Define the custom session class
class SessionWithHeaderRedirection(requests.Session):
    AUTH_HOST = 'urs.earthdata.nasa.gov'

    def __init__(self, username, password):
        super().__init__()
        self.auth = (username, password)

    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers
        url = prepared_request.url

        if 'Authorization' in headers:
            original_parsed = requests.utils.urlparse(response.request.url)
            redirect_parsed = requests.utils.urlparse(url)

            if (original_parsed.hostname != redirect_parsed.hostname) and \
                    redirect_parsed.hostname != self.AUTH_HOST and \
                    original_parsed.hostname != self.AUTH_HOST:
                del headers['Authorization']
        return

# Function to download file
def download_file(session, url, save_dir):
    filename = url[url.rfind('/') + 1:]
    file_path = os.path.join(save_dir, filename)

    try:
        response = session.get(url, stream=True)
        response.raise_for_status()

        with open(file_path, 'wb') as fd:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                fd.write(chunk)
        
        print(f"Downloaded: {file_path}")

    except requests.exceptions.HTTPError as e:
        print(f"HTTP Error for {url}: {e}")
    except Exception as e:
        print(f"Error for {url}: {e}")

doi = '10.5067/EMIT/EMITL1BRAD.001'

# CMR API base url
cmrurl='https://cmr.earthdata.nasa.gov/search/' 

doisearch = cmrurl + 'collections.json?doi=' + doi
concept_id = requests.get(doisearch).json()['feed']['entry'][0]['id']
print(concept_id)



# CMR formatted start and end times
dt_format = '%Y-%m-%dT%H:%M:%SZ'
temporal_str = start_date.strftime(dt_format) + ',' + end_date.strftime(dt_format)
print(temporal_str)

# Search using a Point


point_str = str(lon) +','+ str(lat)

page_num = 1
page_size = 2000 # CMR page size limit

granule_arr = []

while True:
    
     # defining parameters
    cmr_param = {
        "collection_concept_id": concept_id, 
        "page_size": page_size,
        "page_num": page_num,
        "temporal": temporal_str,
        "point":point_str
    }

    granulesearch = cmrurl + 'granules.json'
    response = requests.post(granulesearch, data=cmr_param)
    granules = response.json()['feed']['entry']
       
    if granules:
        for g in granules:
            granule_urls = ''
            granule_poly = ''
                       
            # read cloud cover
            cloud_cover = g['cloud_cover']
    
            # reading bounding geometries
            if 'polygons' in g:
                polygons= g['polygons']
                multipolygons = []
                for poly in polygons:
                    i=iter(poly[0].split (" "))
                    ltln = list(map(" ".join,zip(i,i)))
                    multipolygons.append(Polygon([[float(p.split(" ")[1]), float(p.split(" ")[0])] for p in ltln]))
                granule_poly = MultiPolygon(multipolygons)
            
            # Get https URLs to .nc files and exclude .dmrpp files
            granule_urls = [x['href'] for x in g['links'] if 'https' in x['href'] and '.nc' in x['href'] and '.dmrpp' not in x['href']]
            # Add to list
            granule_arr.append([granule_urls, cloud_cover, granule_poly])
                           
        page_num += 1
    else: 
        break

# creating a pandas dataframe
cmr_results_df = pd.DataFrame(granule_arr, columns=["URL", "cloud_cover", "granule_poly"])
# Drop granules with empty geometry - if any exist
cmr_results_df = cmr_results_df[cmr_results_df['granule_poly'] != '']
# Expand so each row contains a single url 
cmr_results_df = cmr_results_df.explode('URL')
# Name each asset based on filename
cmr_results_df.insert(0,'asset_name', cmr_results_df.URL.str.split('/',n=-1).str.get(-1))

cmr_results_df.to_csv('data/EMIT_extract.csv', index=False)

username = keyring.get_password('emit_service', 'username')
password = keyring.get_password('emit_service', 'password')


# Create a session with the user credentials
session = SessionWithHeaderRedirection(username, password)

# Path to the CSV file and the directory to save files
csv_file_path = "data/EMIT_extract.csv"  # Update this to the path of your CSV file
save_directory = "rads"  # Update this to your desired directory

# Create the directory if it does not exist
os.makedirs(save_directory, exist_ok=True)

# Read the CSV file
df = pd.read_csv(csv_file_path)

# Ensure the URL column exists in the DataFrame
if 'URL' in df.columns:
    for url in df['URL']:
        if 'RAD' in url:
            download_file(session, url, save_directory)
        else:
            print(f"Skipping {url}: does not contain 'RAD'.")
else:
    print("The CSV file does not contain a 'URL' column.")