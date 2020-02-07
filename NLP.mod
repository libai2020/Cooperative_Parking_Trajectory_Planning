param Nfe == 100;
param Nv == 5;
param tf == 40;
param hi = tf / Nfe;
param Nobs := 1;
param BV{i in {1..Nv}, j in {1..6}};
param OV{i in {1..Nobs}, j in {1..Nfe}, k in {1..7}, m in {1..2}};
param Nov{i in {1..Nobs}};
param Area{i in {1..Nobs}};
param amax == 0.5;
param vmax == 2.5;
param wmax == 0.5;
param phymax == 0.7;
param AreaVehicle == 9.106;
var x{i in {1..Nv}, j in {1..Nfe}};
var y{i in {1..Nv}, j in {1..Nfe}};
var xf{i in {1..Nv}, j in {1..Nfe}};
var yf{i in {1..Nv}, j in {1..Nfe}};
var xr{i in {1..Nv}, j in {1..Nfe}};
var yr{i in {1..Nv}, j in {1..Nfe}};
var theta{i in {1..Nv}, j in {1..Nfe}};
var v{i in {1..Nv}, j in {1..Nfe}};
var a{i in {1..Nv}, j in {1..Nfe}};
var phy{i in {1..Nv}, j in {1..Nfe}};
var w{i in {1..Nv}, j in {1..Nfe}};
var egoV{i in {1..Nv}, j in {1..Nfe}, k in {1..4}, m in {1..2}};
minimize cost_function:
sum{i in {1..Nv}, j in {1..Nfe}}(a[i,j]^2); 

s.t. DIFF_dxdt {k in {1..Nv}, i in {2..Nfe}}:
x[k,i] = x[k,i-1] + hi * v[k,i] * cos(theta[k,i]);
s.t. DIFF_dydt {k in {1..Nv}, i in {2..Nfe}}:
y[k,i] = y[k,i-1] + hi * v[k,i] * sin(theta[k,i]);
s.t. DIFF_dvdt {k in {1..Nv}, i in {2..Nfe}}:
v[k,i] = v[k,i-1] + hi * a[k,i];
s.t. DIFF_dthetadt {k in {1..Nv}, i in {2..Nfe}}:
theta[k,i] = theta[k,i-1] + hi * tan(phy[k,i]) * v[k,i] / 2.8;
s.t. DIFF_dphydt {k in {1..Nv}, i in {2..Nfe}}:
phy[k,i] = phy[k,i-1] + hi * w[k,i];

s.t. RELATIONSHIP_AX {k in {1..Nv}, i in {1..Nfe}}:
egoV[k,i,1,1] = x[k,i] + 3.76 * cos(theta[k,i]) - 0.971 * sin(theta[k,i]);
s.t. RELATIONSHIP_BX {k in {1..Nv}, i in {1..Nfe}}:
egoV[k,i,2,1] = x[k,i] + 3.76 * cos(theta[k,i]) + 0.971 * sin(theta[k,i]);
s.t. RELATIONSHIP_CX {k in {1..Nv}, i in {1..Nfe}}:
egoV[k,i,3,1] = x[k,i] - 0.929 * cos(theta[k,i]) + 0.971 * sin(theta[k,i]);
s.t. RELATIONSHIP_DX {k in {1..Nv}, i in {1..Nfe}}:
egoV[k,i,4,1] = x[k,i] - 0.929 * cos(theta[k,i]) - 0.971 * sin(theta[k,i]);
s.t. RELATIONSHIP_AY {k in {1..Nv}, i in {1..Nfe}}:
egoV[k,i,1,2] = y[k,i] + 3.76 * sin(theta[k,i]) + 0.971 * cos(theta[k,i]);
s.t. RELATIONSHIP_BY {k in {1..Nv}, i in {1..Nfe}}:
egoV[k,i,2,2] = y[k,i] + 3.76 * sin(theta[k,i]) - 0.971 * cos(theta[k,i]);
s.t. RELATIONSHIP_CY {k in {1..Nv}, i in {1..Nfe}}:
egoV[k,i,3,2] = y[k,i] - 0.929 * sin(theta[k,i]) - 0.971 * cos(theta[k,i]);
s.t. RELATIONSHIP_DY {k in {1..Nv}, i in {1..Nfe}}:
egoV[k,i,4,2] = y[k,i] - 0.929 * sin(theta[k,i]) + 0.971 * cos(theta[k,i]);

s.t. RELATIONSHIP_XF {k in {1..Nv}, i in {1..Nfe}}:
xf[k,i] = x[k,i] + 2.5877 * cos(theta[k,i]);
s.t. RELATIONSHIP_YF {k in {1..Nv}, i in {1..Nfe}}:
yf[k,i] = y[k,i] + 2.5877 * sin(theta[k,i]);
s.t. RELATIONSHIP_XR {k in {1..Nv}, i in {1..Nfe}}:
xr[k,i] = x[k,i] + 0.2432 * cos(theta[k,i]);
s.t. RELATIONSHIP_YR {k in {1..Nv}, i in {1..Nfe}}:
yr[k,i] = y[k,i] + 0.2432 * sin(theta[k,i]);

s.t. Bonds_phy {k in {1..Nv}, i in {1..Nfe}}:
-phymax <= phy[k,i] <= phymax;
s.t. Bonds_a {k in {1..Nv}, i in {1..Nfe}}:
-amax <= a[k,i] <= amax;
s.t. Bonds_v {k in {1..Nv}, i in {1..Nfe}}:
-vmax <= v[k,i] <= vmax;
s.t. Bonds_w {k in {1..Nv}, i in {1..Nfe}}:
-wmax <= w[k,i] <= wmax;
s.t. Bonds_x {k in {1..Nv}, i in {1..Nfe}}:
-20 <= x[k,i] <= 20;
s.t. Bonds_y {k in {1..Nv}, i in {1..Nfe}}:
-20 <= y[k,i] <= 20;

s.t. ObsVertexOutOfABCD {m in {1..Nv}, i in {1..Nfe}, j in {1..Nobs}, k in {1..Nov[j]}}:
0.5 * abs(OV[j,i,k,1] * egoV[m,i,1,2] + egoV[m,i,1,1] * egoV[m,i,2,2] + egoV[m,i,2,1] * OV[j,i,k,2] - OV[j,i,k,1] * egoV[m,i,2,2] - egoV[m,i,1,1] * OV[j,i,k,2] - egoV[m,i,2,1] * egoV[m,i,1,2]) + 
0.5 * abs(OV[j,i,k,1] * egoV[m,i,3,2] + egoV[m,i,3,1] * egoV[m,i,2,2] + egoV[m,i,2,1] * OV[j,i,k,2] - OV[j,i,k,1] * egoV[m,i,2,2] - egoV[m,i,3,1] * OV[j,i,k,2] - egoV[m,i,2,1] * egoV[m,i,3,2]) + 
0.5 * abs(OV[j,i,k,1] * egoV[m,i,3,2] + egoV[m,i,3,1] * egoV[m,i,4,2] + egoV[m,i,4,1] * OV[j,i,k,2] - OV[j,i,k,1] * egoV[m,i,4,2] - egoV[m,i,3,1] * OV[j,i,k,2] - egoV[m,i,4,1] * egoV[m,i,3,2]) + 
0.5 * abs(OV[j,i,k,1] * egoV[m,i,1,2] + egoV[m,i,1,1] * egoV[m,i,4,2] + egoV[m,i,4,1] * OV[j,i,k,2] - OV[j,i,k,1] * egoV[m,i,4,2] - egoV[m,i,1,1] * OV[j,i,k,2] - egoV[m,i,4,1] * egoV[m,i,1,2]) >= AreaVehicle + 0.1;

s.t. CarVertexOutOfObstacle {m in {1..Nv}, i in {1..Nfe}, j in {1..Nobs}, k in {1..4}}:
0.5 * abs(egoV[m,i,k,1] * OV[j,i,1,2] + OV[j,i,1,1] * OV[j,i,2,2] + OV[j,i,2,1] * egoV[m,i,k,2] - egoV[m,i,k,1] * OV[j,i,2,2] - OV[j,i,1,1] * egoV[m,i,k,2] - OV[j,i,2,1] * OV[j,i,1,2]) + 
0.5 * abs(egoV[m,i,k,1] * OV[j,i,3,2] + OV[j,i,3,1] * OV[j,i,2,2] + OV[j,i,2,1] * egoV[m,i,k,2] - egoV[m,i,k,1] * OV[j,i,2,2] - OV[j,i,3,1] * egoV[m,i,k,2] - OV[j,i,2,1] * OV[j,i,3,2]) + 
0.5 * abs(egoV[m,i,k,1] * OV[j,i,3,2] + OV[j,i,3,1] * OV[j,i,4,2] + OV[j,i,4,1] * egoV[m,i,k,2] - egoV[m,i,k,1] * OV[j,i,4,2] - OV[j,i,3,1] * egoV[m,i,k,2] - OV[j,i,4,1] * OV[j,i,3,2]) + 
0.5 * abs(egoV[m,i,k,1] * OV[j,i,1,2] + OV[j,i,1,1] * OV[j,i,4,2] + OV[j,i,4,1] * egoV[m,i,k,2] - egoV[m,i,k,1] * OV[j,i,4,2] - OV[j,i,1,1] * egoV[m,i,k,2] - OV[j,i,4,1] * OV[j,i,1,2]) >= Area[j] + 0.1;

s.t. VehicleItoJff {i in {1..(Nv-1)}, j in {(i+1)..Nv}, kk in {1..Nfe}}:
(xf[i,kk] - xf[j,kk])^2 + (yf[i,kk] - yf[j,kk])^2 >= 9.2680;
s.t. VehicleItoJrr {i in {1..(Nv-1)}, j in {(i+1)..Nv}, kk in {1..Nfe}}:
(xr[i,kk] - xr[j,kk])^2 + (yr[i,kk] - yr[j,kk])^2 >= 9.2680;
s.t. VehicleItoJfr {i in {1..(Nv-1)}, j in {(i+1)..Nv}, kk in {1..Nfe}}:
(xf[i,kk] - xr[j,kk])^2 + (yf[i,kk] - yr[j,kk])^2 >= 9.2680;
s.t. VehicleItoJrf {i in {1..(Nv-1)}, j in {(i+1)..Nv}, kk in {1..Nfe}}:
(xr[i,kk] - xf[j,kk])^2 + (yr[i,kk] - yf[j,kk])^2 >= 9.2680;

s.t. EQ_init_x {m in {1..Nv}}:
x[m,1] = BV[m,1];
s.t. EQ_init_y {m in {1..Nv}}:
y[m,1] = BV[m,2];
s.t. EQ_init_theta {m in {1..Nv}}:
theta[m,1] = BV[m,3];
s.t. EQ_init_v {m in {1..Nv}}:
v[m,1] = 0;
s.t. EQ_init_a {m in {1..Nv}}:
a[m,1] = 0;
s.t. EQ_init_phy {m in {1..Nv}}:
phy[m,1] = 0;
s.t. EQ_init_w {m in {1..Nv}}:
w[m,1] = 0;
s.t. EQ_end_x {m in {1..Nv}}:
x[m,Nfe] = BV[m,4];
s.t. EQ_end_y {m in {1..Nv}}:
y[m,Nfe] = BV[m,5];
s.t. EQ_end_theta1 {m in {1..Nv}}:
sin(theta[m,Nfe]) = sin(BV[m,6]);
s.t. EQ_end_theta2 {m in {1..Nv}}:
cos(theta[m,Nfe]) = cos(BV[m,6]);
s.t. EQ_end_v {m in {1..Nv}}:
v[m,Nfe] = 0;
s.t. EQ_end_a {m in {1..Nv}}:
a[m,Nfe] = 0;
s.t. EQ_end_phy {m in {1..Nv}}:
phy[m,Nfe] = 0;
s.t. EQ_end_w {m in {1..Nv}}:
w[m,Nfe] = 0;

data;
param: BV := include BV;
param: OV := include OV;
param: Nov := include Nov;
param: Area := include Area;