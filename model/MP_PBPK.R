mPBPK_RES <- '
$PARAM
// ---------- Physiological basics ----------
BW   = 0.02        // kg
QCC  = 16.5        // L/h/kg^0.75

// ---------- Organ blood flow fractions (of QC) ----------
f_QGI = 0.141
f_QL  = 0.161
f_QK  = 0.091
f_QS  = 0.011
f_QLU = 1

// ---------- Volume fractions (of BW; units: L since BW in kg) ----------
v_VGI = 0.042
v_VL  = 0.055
v_VK  = 0.017
v_VS  = 0.005
v_VB  = 0.049
v_VLU = 0.007

// ---------- Vascular fractions (within each tissue) ----------
BVGI = 0.03
BVL  = 0.31
BVK  = 0.24
BVS  = 0.17
BVLU = 0.50
BVR  = 0.04      // vascular fraction in Rest

// ---------- Permeability-surface (as flow fractions) ----------
PAGI = 0.001
PAL  = 0.001
PAK  = 0.001
PAS  = 0.001
PALU = 0.001
PAR  = 1e-4      

// ---------- Elimination & GI loop ----------
KbileC  = 0.0022    // L/h/kg^0.75, Biliary clearance; Lin et al., 2016
KurineC = 0.00032  // L/h/kg^0.75, Urinary clearance; Lin et al., 2016
Kfeces = 0.53     // h^-1 from lumen

// ---------- RES kinetics (capacity-limited uptake; first-order release) ----------
// Liver RES
KLRESmax     = 0.4      // h^-1 (maximum uptake scaling)
KLRESrelease = 0.0075   // h^-1 release from RES to liver tissue/blood
ALREScap     = 100    // mg/L * L = mg capacity scale (cap * VL)

// Spleen RES
KSRESmax     = 1
KSRESrelease = 0.003
ASREScap     = 100

// Kidney RES
KKRESmax     = 0.1
KKRESrelease = 0.01
AKREScap     = 15

// Lung RES
KLURESmax     = 0.1
KLURESrelease = 0.005
ALUREScap     = 15

// Rest of Body RES
KRRESmax     = 0.2     // h^-1  max uptake scale in rest
KRRESrelease = 0.002   // h^-1  slow release from rest RES
ARREScap     = 15    // mg capacity scale (cap * VR)


// ---------- Partition coefficients (tissue:blood, unitless) ----------
KpL  = 0.01   // liver
KpK  = 0.15   // kidney
KpS  = 0.15   // spleen
KpGI = 0.15   // GI tissue (not lumen)
KpLU = 0.15 
KpR  = 0.15   // rest-of-body

// ---------- Gut-wall barrier kinetics ----------
Ktrans = 5e-1   // h^-1  lumen -> gut-wall
Kgw2b  = 3e-4   // h^-1  gut-wall -> GI blood
Kshed  = 5e-1   // h^-1  gut-wall -> lumen (shedding)

// #+ Albumin-binding parameters
KD                 = 206.15    
N                  = 1         
Calb               = 0.035     
Mwalb              = 66500     
MwPS               = 191600    

$INIT
A_dose      = 0
A_blood     = 0
A_liver_b   = 0
A_liver_t   = 0
A_LRES      = 0   // Liver RES

A_kidney_b  = 0
A_kidney_t  = 0
A_KRES      = 0   // Kidney RES

A_spleen_b  = 0
A_spleen_t  = 0
A_SRES      = 0   // Spleen RES

A_GI_b      = 0
A_GI_t      = 0
A_lumen     = 0
A_urine     = 0
A_bile      = 0
A_feces     = 0

A_lung_b    = 0
A_lung_t    = 0
A_LURES     = 0

A_rest_b    = 0
A_rest_t    = 0
      

$MAIN
// ---- Cardiac output (L/h) ----
double QC  = QCC * pow(BW, 0.75);
double QGI = QC * f_QGI;
double QL  = QC * f_QL;
double QK  = QC * f_QK;
double QS  = QC * f_QS;
double QLU = QC * f_QLU;  
double QR  = QC * (1 - (f_QGI + f_QL + f_QK + f_QS));
double QBal = QC - (QGI + QL + QK + QS + QR);

// ---- Volumes (L) ----
double VGI = BW * v_VGI;
double VL  = BW * v_VL;
double VK  = BW * v_VK;
double VS  = BW * v_VS;
double VB  = BW * v_VB;
double VLU = BW*v_VLU;
double VR  = BW - (VGI + VL + VK + VS + VB + VLU);
double VBal = BW - (VGI + VL + VK + VS + VB + VLU + VR);

// ---- Sub-volumes (blood/tissue splits) ----
double VGIb = VGI * BVGI; double VGIt = VGI - VGIb;
double VLb  = VL  * BVL;  double VLt  = VL  - VLb;
double VKb  = VK  * BVK;  double VKt  = VK  - VKb;
double VSb  = VS  * BVS;  double VSt  = VS  - VSb;
double VLUb = VLU * BVLU; double VLUt = VLU - VLUb;
double VRb = VR * BVR; double VRt = VR - VRb;

// ---- PS "flows" (L/h) ----
double PAGI_flow = PAGI * QGI;
double PAL_flow  = PAL  * QL;
double PAK_flow  = PAK  * QK;
double PAS_flow  = PAS  * QS;
double PALU_flow = PALU * QLU;
double PAR_flow = PAR * QR;   // “PS flow” for Rest

// #+ Biliary excretion 
double Kbile       = KbileC*pow(BW, 0.75);
double Kurine      = KurineC*pow(BW, 0.75);

// #+ Maximum protein-binding capacity (mg/L)
double Bmax     = Calb*N*MwPS/Mwalb*1000;

// #+ Unbound fraction in blood (-)
double fu      = 1; //(KD)/(Bmax+KD);


$ODE
// ============================================================
// 1) Dose input: A_dose -> A_blood (first-order)
// ============================================================
dxdt_A_dose = 0;

// ============================================================
// 2) Concentrations (mg/L)
// ============================================================
double CA   = A_blood   / VB;
double CLb  = A_liver_b / VLb;
double CLt  = A_liver_t / VLt;
double CKb  = A_kidney_b/ VKb;
double CKt  = A_kidney_t/ VKt;
double CSb  = A_spleen_b/ VSb;
double CSt  = A_spleen_t/ VSt;
double CGIb = A_GI_b    / VGIb;
double CGIt = A_GI_t    / VGIt;



double CLUb = A_lung_b  / VLUb;   
double CLUt = A_lung_t  / VLUt;   

double CRb = A_rest_b / VRb;
double CRt = A_rest_t / VRt;

// ============================================================
// 3) Blood (arterial/venous lumped as A_blood here):
//    Distributes to organs and receives venous return
// ============================================================
dxdt_A_blood = QL*(CLb - CA)
             + QK*(CKb - CA)
             + QGI*(CGIb - CA)
             + QS*(CSb - CA)
             + QLU*(CLUb - CA)
             + QR*(CRb - CA);

// ============================================================
// 4) Liver with RES (uptake from liver blood, release to liver tissue)
//    Uptake capacity: ALREScap * VL
// ============================================================
double KLRESUP = KLRESmax * (1.0 - A_LRES/(ALREScap * VL));
if (KLRESUP < 0.0) KLRESUP = 0.0;

double J_L = PAL_flow*(CLb - CLt/KpL)*fu;

dxdt_A_liver_b = QL*(CA - CLb) - J_L - KLRESUP*A_liver_b;
                
dxdt_A_liver_t = J_L  - Kbile*CLt
                 + KLRESrelease*A_LRES;                // RES release back 
                
dxdt_A_LRES    = KLRESUP*A_liver_b - KLRESrelease*A_LRES;  // RES balance (internal transfer)

dxdt_A_bile    = Kbile*CLt;   // biliary excretion from liver tissue

// ============================================================
// 5) Kidney with RES (uptake from kidney tissue, release to kidney tissue)
//    Uptake capacity: AKREScap * VK
// ============================================================
double KKRESUP = KKRESmax * (1.0 - A_KRES/(AKREScap * VK));
if (KKRESUP < 0.0) KKRESUP = 0.0;

double J_K = PAK_flow*(CKb - CKt/KpK)*fu;

dxdt_A_kidney_b = QK*(CA - CKb) - J_K - Kurine*CKb*fu;

dxdt_A_kidney_t = J_K
                  - KKRESUP*A_kidney_t
                  + KKRESrelease*A_KRES;             // RES release to tissue

dxdt_A_KRES     = KKRESUP*A_kidney_t - KKRESrelease*A_KRES;

dxdt_A_urine    = Kurine*CKb*fu;   // urinary excretion

// ============================================================
// 6) Spleen with RES (uptake from spleen tissue, release to spleen tissue)
//    Uptake capacity: ASREScap * VS
// ============================================================
double KSRESUP = KSRESmax * (1.0 - A_SRES/(ASREScap * VS));
if (KSRESUP < 0.0) KSRESUP = 0.0;

double J_S = PAS_flow*(CSb - CSt/KpS)*fu;

dxdt_A_spleen_b = QS*(CA - CSb) - J_S;

dxdt_A_spleen_t = J_S
                  - KSRESUP*A_spleen_t
                  + KSRESrelease*A_SRES;            // RES release to tissue

dxdt_A_SRES     = KSRESUP*A_spleen_t - KSRESrelease*A_SRES;

// ============================================================
// 7) GI tract with enterohepatic recirculation
// ============================================================
double Ka_eff = (Ktrans * Kgw2b) / (Kgw2b + Kshed);
double J_GI = PAGI_flow*(CGIb - CGIt/KpGI)*fu;

dxdt_A_GI_b = QGI*(CA - CGIb)
              - J_GI
              + Ka_eff * A_lumen;                       // reabsorption

dxdt_A_GI_t = J_GI;

dxdt_A_lumen = Kbile*CLt - (Kfeces + Ka_eff)*A_lumen;

dxdt_A_feces = Kfeces*A_lumen;


// ============================================================
// 7) Lung 
// ============================================================
double KLURESUP = KLURESmax * (1.0 - A_LURES/(ALUREScap * VLU));
if (KLURESUP < 0.0) KLURESUP = 0.0;

double J_LU = PALU_flow*(CLUb - CLUt/KpLU)*fu;

dxdt_A_lung_b = QLU*(CA - CLUb) - J_LU;
dxdt_A_lung_t = J_LU
                - KLURESUP*A_lung_t
                + KLURESrelease*A_LURES;
                
dxdt_A_LURES = KLURESUP*A_lung_t - KLURESrelease*A_LURES;

// ============================================================
// 9) Rest-of-body: flow-limited tissue (no PS term)
// ============================================================

double J_R = PAR_flow * (CRb - CRt / KpR) * fu;

dxdt_A_rest_b = QR*(CA - CRb) - J_R;
dxdt_A_rest_t = J_R;

// ============================================================
// 10) Mass balance (system + excreted vs. input)
// ============================================================
double A_body = A_blood 
              + A_liver_b + A_liver_t + A_LRES
              + A_kidney_b + A_kidney_t + A_KRES
              + A_spleen_b + A_spleen_t + A_SRES
              + A_GI_b + A_GI_t + A_lumen
              + A_lung_b + A_lung_t + A_LURES
              + A_rest_b + A_rest_t;            

double A_excreted = A_urine + A_feces;

double A_total = A_body + A_excreted;

double mass_balance = A_dose - A_total;

double total_liver  = A_liver_b + A_liver_t + A_LRES;
double total_kidney = A_kidney_b + A_kidney_t + A_KRES;
double total_spleen = A_spleen_b + A_spleen_t + A_SRES;
double total_lung   = A_lung_b + A_lung_t + A_LURES;
double total_GI     = A_GI_b + A_GI_t + A_lumen;
double total_rest = A_rest_b + A_rest_t;

$CAPTURE A_body A_excreted A_total mass_balance
         total_liver total_kidney total_spleen total_GI
         total_lung total_rest
'
