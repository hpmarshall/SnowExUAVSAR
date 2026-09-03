# SnowExUAVSAR
this repo summarizes the UAVSAR data collected for SnowEx, and starts by building a database for Mores Creek Summit area

## Tools to use
Use asf_search for accessing UAVSAR data
Use uavsar_pytools for reading the binary UAVSAR data

## Area of interest
Mores Creek Summit, MCS_domain.kml

## First task
Start by querying the UAVSAR data on ASF using asf_search, and find all UAVSAR flights during the periods January 1, 2020-April 1, 2020 and January 1, 2021-April 1, 2021, in the Western U.S.

Make a table of all flights with location, state, date, and bearing.

## Plan Mode
Always use plan mode and create a plan for tasks, and ask for permission before executing.

## Code style
Keep code as concise as possible, and as simple to understand.  Comment every line. 

## References
Primary reference for SnowEx UAVSAR data: https://nsidc.org/sites/default/files/documents/technical-reference/snex_uavsar-v001-techref.pdf

Primary references for SnowEx LiDAR:
https://nsidc.org/data/snex20_qsi_sd_3m/versions/1
https://nsidc.org/data/snex23_lidar/versions/1
https://nsidc.org/sites/default/files/documents/user-guide/snex21_ps_dsm-v001-userguide.pdf
https://nsidc.org/sites/default/files/documents/user-guide/snex20_gm_swe_sd.pdf

## Extended search 
SnowEx 2017 (Grand Mesa, Telluride)  and flights in 2014-2016 in mountainous areas

## Task 2
Create kml of all lidar acquisitions associated with SnowEx.  These should come from the references for SnowEx LiDAR above, and the Airborne Snow Observatory.  We are interested in any airborne lidar flown during the months January-April that overlaps with the UAVSAR flights from 2014-2021.
