# StormWISE model Cost minimization with multiple benefits as constraints 
# Adapted for Wingohocking with slopes obtained from SWMM runs
# Replace landuse category with land ownership category   
# Drainage zone is subcatchment based on SWMM analysis
set I;	# drainage zones
set J;  # ownership category
set K;	# bmp/lid categories 
set T;	# benefit categories
set D;  # deployment categories

param Bmin{T}; # default 0.0;		# lower bounds for benefits in converted units
param convert{T}; # conversion factors for benefits
param deploy{K,D};  # fraction deployable on deployment category d
param pctown{I,J};  # percent of land in each ownership category
# treatment fraction calculations:
param impfr{I};  # impervious fraction zone i
param f{i in I,j in J,k in K} =  deploy[k,"pervious"]*(1-impfr[i]) + deploy[k,"impervious"]*impfr[i];	
param area{I};	# area in zone i 
param h{I,K,T};   # from SWMM - slopes of volume reduction vs number of LID
param g{K};      # greened acre per one unit of GI type k 
param a{J};       # AKRF cost model intercept
param b{J};       # AKRF cost model slope
param cost{j in J, k in K} = 10**(a[j] - b[j]*log10(g[k])); # GI practice cost per greened Acreparam h{I,K,T};      
# amount of benefit T per one unit of GI type k in zone i (from SWMM for runoff):
param s{i in I,k in K,t in T} = h[i,k,t]/g[k];	# calculated benefit slopes
param u{i in I,j in J,k in K} = f[i,j,k]*area[i];	# calculated upper bnds on greened acres

var x{i in I,j in J,k in K} >= 0 <= u[i,j,k];	# number of greened acres - Decision Variables

minimize investment: sum{i in I,j in J,k in K} cost[j,k]*x[i,j,k];

subject to benefits{t in T}: sum{i in I, j in J, k in K} s[i,k,t]*x[i,j,k] >= Bmin[t]/convert[t];
