RUNSPEC
TITLE
CO2 storage simulation on Smeaheia
--
-- Pore volume multipliers in North and South
--
-- Run details
--  Dissolved gas: NO
--  Pc: NO
--  Rel.perm: Linear, Sor=0.15, Sgr=0.05
 
 
 
NOECHO
 
DIMENS
-- Grid Dimensions
-- NX	NY	NZ
-- ------  --------  ------
   33  60  27  /
 
-- Active Phases Present
OIL
GAS
--DISGAS
CO2STORE

METRIC
-- Unit Convention
 
--BIGMODEL
-- DIFFUSE
-- Enables Molecular Diffusion
 
PARALLEL
 4 DISTRIBUTED /
 
--MEMORY
-- 1000 /
 
TABDIMS
-- Table Of Dimensions
-- NTSFUN  NTPVT  NSSFUN  NPPVT   NTFIP
-- ------  -----  ------ ----- -----
     2       1      30    30    8   /
 
-- NTSFUN: No. of saturation tables entered.
-- NTPVT : No. of PVT tables entered (in the PROPS section).
-- NSSFUN: Max. no. of saturation node in each saturation table, ie.,
--         Max. no. of data points in each table.
 
WELLDIMS
-- Well Dimension Data
-- NWMAXZ NCWMAX NGMAXZ NWGMAX
-- ------ ------ ------ ------
     145      30    3    4     /
 
-- NWMAXZ: Max. no. of wells in the models.
-- NCWMAX: Max. no. of conncections per well (i.e., no. of perforations).
-- NGMAXZ: Max. no. of groups in the model.
-- NWGMAX: Max. no. of wells in any group.
 
AQUDIMS
2 10000/


START
-- Spesifies a Start Date
-- DAY    MONTH   YEAR
-- ----   -----   ----
    1     JAN    1991 /
 
FAULTDIM
300 /
 
NSTACK
-- Stack Size For Linear Solver
 50 /
 
UNIFOUT
-- Restart And Summary Files Are To Be Unified
 
UNIFIN
-- Restart From A Unified Restart File
 
--NOSIM

--GRIDOPTS
-- YES
--/
GRID

GRIDFILE
 0 0 /

INCLUDE
 '../GRID_BOX_${LEVEL}.INC'
/

INCLUDE
 'smeaheia_${ind}.perm' 
/

--POROSITY
--INCLUDE
--'../include/utsira_200_10layers_poro.dat' /
 
--PORO
-- 418338*0.35 /
 
--EQUALS
-- 'PORO' 0.01 1 127 1 183   6  6 /
-- 'PORO' 0.01 1 127 1 183  10 10 /
-- 'PORO' 0.01 1 127 1 183  14 14 /
-- 'PORO' 0.01 1 127 1 183  18 18 /
--/
 
-- Removing eastern block to get rid of pvt problem at low pressure
--EQUALS
--  ACTNUM   0   66 109  75 238 1 27 /
--/
 
BOX
${PORO_I} ${PORO_I} 1 60 1 27 /
--33 33 1 60 1 27 /
--65 65 1 60 1 27 /
--97 97 1 60 1 27 /

PORO 
 1620*10000/

ENDBOX


PINCH
0.1 /
 
MAPAXES
0.0 100.0 0.0 0.0 0.0 0.0 /
 
INCLUDE
 vette_${ind}.aqucon
/

AQUNUM
-- id  I J K Area Length Poro K Depth Initial Pr PVTNUM SATNUM 
1 1 60 1 40 500 0.1 ${K1} 1* 1* 1 2 /
1 2 60 1 40 500 1e12 1e3 1* 1* 1 1 /
2 1 58 1 4000 0.1 1e12 ${K2} 1100 110 1 1 /
/

INCLUDE
 vette_${ind}.multx /
 
--INIT
-- ===============================================================================
EDIT
 
-- Pore volume multipliers along western and northern edge of model
--INCLUDE
-- '../include/utsira_200_10_PVMULT_NW.dat' /
 
-----------------------------------
--===========================ls ================================================
PROPS
-- ===============================================================================
 

--------------------------------------------------------------------------------
-- Start include file ../include/UTS_S_DISGAS_rc4-5_ext.ECL
--------------------------------------------------------------------------------


-- =======================
-- Pressure Dependent Data
-- =======================
 
ROCK
-- Rock Compressibility
-- Ref. pressure     Compressibility
-- -------------     ---------------
       10.0              4E-05   /
 


INCLUDE
 vette_${ind}.satfun /
 
-- ===============================================================================
REGIONS
-- ==============================================================================
FIPNUM
${FIPNUM_NUM}*1 /
--53460*1 /
--105300*1 /
--157140*1 /


BOX
1 2 60 60 1 1 /
FIPNUM
 2*2/

ENDBOX

BOX
1 1 58 58 1 1 /
FIPNUM
 1*3/

ENDBOX

-- ===============================================================================
SOLUTION
-- ===============================================================================
 
EQUIL
-- Equilibration Data Spesification
-- Datum  Pi@Datum  WOC   Pc@WOC   GOC   Pc@GOC  Rs  Rv  Accuracy
-- -----  --------  ---   ------  -----  ------  --  --  --------
    800     81.0   5050.0    0.0    400.0    0.0    1   0          /
 
 
RPTSOL
-- Report Switches For SOLUTION Data
--  p   So  Sw  Sg  Ps  xm  rst   vm    ml  rs
--  1    2   3   4   5   6   7     8    9   10
-- ---  --- --- --- --- --- ---   ---  ---  --
--    1    1   0   1   1   0   3     1    4   1   /
 FIP=2 RESTART=2 /
 
 
RSVD
-- Variation Of Solution GOR With Depth
-- Depth     Rs
-- -----  --------
    800   0.00000
   4150   0.00000 /
 
--RPTRST
--DEN /
 
-- ===============================================================================
SUMMARY
-- ===============================================================================
ECHO
 
 
FPR
 
FOPR
 
FOPT
 
FGPR
 
FGPT
 
FOIR
 
FOIT
 
FGIR
 
FGIT
 
FOIP
 
FGIP
 
WBHP
/
 
TCPU
 
ELAPSED
 
-- Injector
 
FVIR
 
FVIT
 
 
ANQP
/

FNQR 
/

ANQR
/

ANQT
/

ROIP
/

RGIP
/
RGIPL
/
RGIPG
/

FGIP
/

FOIP
/

RPTONLY
RUNSUM
EXCEL
 
-- ===============================================================================
SCHEDULE
-- ===============================================================================
 
 
MESSAGES
-- Resets Message Print and Stop Limits
-- Messages   Comments  Warnings  Problems  Error    Bug
-- --------   --------  --------  --------  -----    ---
    1000  1000 1000000 1000  100 100 10000000 10000000 10000000 1000000 /
 
 
RPTSCHED
'FIP=3' 'CPU=1' 'WELSPECS' /
--'FIP=3' 'WELLS=5' 'SUMMARY=1' 'CPU=1' 'WELSPECS' 'NEWTON=1' /
-- 'NOTHING' /
-- 'RESTART=2' 'RS' 'SGAS' 'SOIL' 'PRESSURE' /
-- Report Switches For Simulation Results
--                        Restart files at
--  P  So  Sw  Sg  Rs  Rv every report time  FIP  Well VFPPROD Sum. CPU    NEWTON
-- --- --  --  --  --  -- -----------------  ---  ---- ------- ---  --- -- ------
--    1   1   0   1   1   0        2            1     0     0    0    2   0    1  /
 
RPTRST
 BASIC=0 /
--BASIC=3 FREQ=2 /
--3 4* 2  /
-- BASIC=2 PRES FIP /
 
--TUNING
--1.0 5.0 0.1 0.15 3.0 0.3 0.4 1.5/
--1.0 0.2 5.0E-7 0.005 10.0 10.0 5.0E-6 0.05 0.001 0.5 /
--12 1 40 1 50 8 /
 
DRSDT
0 /
 
 
WELSPECS
 'SDL1' 'INJ'    ${WELL_I}  20  1*  'GAS'  0.2 /
-- 'SDL1' 'INJ'    10   20  1*  'GAS'  0.2 /
-- 'SDL1' 'INJ'    14   20  1*  'GAS'  0.2 /
--  'SDL1' 'INJ'    15   20  1*  'GAS'  0.2 /
/
 
COMPDAT
 'SDL1' ${WELL_I}   20   16  19  'OPEN' 0  1*  0.2 /
-- 'SDL1' 10   20   16  19  'OPEN' 0  1*  0.2 /
--  'SDL1' 14   20   16  19  'OPEN' 0  1*  0.2 /
--    'SDL1' 15   20   16  19  'OPEN' 0  1*  0.2 /
/
 
 
WCONINJE
 'SDL1' 'GAS'    'OPEN'    'RATE'      2202000     1*     300   /
/

--WELSPECS
-- 'TROLL' 'PROD'    1 58   100  'OIL'  0.2 /
--/
 
--COMPDAT
-- 'TROLL' 1 58   1  1  'OPEN' 0  1*  0.2 /
--/
 
 
--WCONPROD
-- 'TROLL' 'OPEN' 'ORAT' 1e3 4* 100 /  
--/

DATES
1 'FEB' 1991 /
1 'JLY' 1991 /
1 'JAN' 1992 /
1 'JAN' 1993 /
1 'JAN' 1994 /
1 'JAN' 1995 /
1 'JAN' 1996 /
1 'JAN' 1997 /
/

DATES
1 'JAN' 1999 /
1 'JAN' 1999 /
1 'JAN' 2000 /
1 'JAN' 2001 /
1 'JAN' 2002 /
1 'JAN' 2003 /
1 'JAN' 2004 /
1 'JAN' 2005 /
1 'JAN' 2006 /
1 'JAN' 2007 /
1 'JAN' 2008 /
1 'JAN' 2009 /
1 'JAN' 2010 /
1 'JAN' 2011 /
1 'JAN' 2012 /
1 'JAN' 2013 /
1 'JAN' 2014 /
1 'JAN' 2015 /
1 'JAN' 2016 /
1 'JAN' 2017 /
1 'JAN' 2018 /
1 'JAN' 2019 /
1 'JAN' 2020 /
1 'JAN' 2021 /
1 'JAN' 2022 /
1 'JAN' 2023 /
1 'JAN' 2024 /
1 'JAN' 2025 /
1 'JAN' 2026 /
1 'JAN' 2027 /
1 'JAN' 2028 /
1 'JAN' 2029 /
1 'JAN' 2030 /
1 'JAN' 2031 /
1 'JAN' 2032 /
1 'JAN' 2033 /
1 'JAN' 2034 /
1 'JAN' 2035 /
1 'JAN' 2036 /
1 'JAN' 2037 /
1 'JAN' 2038 /
1 'JAN' 2039 /
1 'JAN' 2040 /
1 'JAN' 2041 /
1 'JAN' 2042 /
1 'JAN' 2043 /
1 'JAN' 2044 /
1 'JAN' 2045 /
1 'JAN' 2046 /
1 'JAN' 2047 /
1 'JAN' 2048 /
1 'JAN' 2049 /
1 'JAN' 2050 /
/
