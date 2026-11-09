"""
Corrected RSM Table 3 reconstruction and hard-constraint optimisation
for manuscript CSCEE-D-26-00383.

Key correction:
- FWS31 = 1S:3FW
- FWS12 = 2S:1FW
- FWS13 = 3S:1FW
- all three microbial models must satisfy the selected FIB threshold simultaneously
- candidate settings are restricted to the Box-Behnken convex design hull

Run from the repository root:
    python python/corrected_rsm_table3_hard_constraints.py
"""
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
data_path = ROOT / "data" / "input" / "trial_data.csv"
if not data_path.exists():
    raise FileNotFoundError(f"Could not find {data_path}")

df = pd.read_csv(data_path)
recipe_map = {"FWS31":"1S:3FW", "FWS12":"2S:1FW", "FWS13":"3S:1FW"}
df["recipe"] = df["Recipe"].map(recipe_map)

for col in ["Coliforms","Ecoli","Enterococci"]:
    df[col+"_t"] = np.log10(df[col] + 1.0)
    df.loc[df[col+"_t"] < 1.0, col+"_t"] = 1.0

def xrow(temp, olr, rt, recipe):
    i2 = float(recipe == "2S:1FW")
    i3 = float(recipe == "3S:1FW")
    return np.array([
        1,temp,olr,rt,i2,i3,temp**2,olr**2,rt**2,
        temp*olr,temp*rt,olr*rt,temp*i2,temp*i3,
        olr*i2,olr*i3,rt*i2,rt*i3
    ], float)

X = np.vstack([xrow(r.Temp,r.OL,r.Time,r.recipe) for r in df.itertuples()])
responses = {
    "Coliforms":"Coliforms_t",
    "E. coli":"Ecoli_t",
    "Enterococci":"Enterococci_t",
    "Methane":"Methane",
}
models = {name: sm.OLS(df[col].to_numpy(float), X).fit() for name,col in responses.items()}
betas = {name: np.asarray(model.params,float) for name,model in models.items()}

def predict(temp,olr,rt,recipe):
    xx=xrow(temp,olr,rt,recipe)
    return {name:float(xx@beta) for name,beta in betas.items()}

def hull_margin(temp,olr,rt):
    tc=(temp-37)/18
    oc=(olr-2)/1.5
    rc=(rt-6)/3
    return 2-(abs(tc)+abs(oc)+abs(rc))

def optimise(recipe, threshold, objective):
    cons=[{"type":"ineq","fun":lambda z:hull_margin(*z)}]
    for response in ["Coliforms","E. coli","Enterococci"]:
        cons.append({"type":"ineq",
                     "fun":lambda z,r=response:threshold-predict(*z,recipe)[r]})
    fun=(lambda z:z[0]) if objective=="min_temp" else (lambda z:-predict(*z,recipe)["Methane"])
    starts=[[48.5,.5,7],[50.5,2.3,9],[51.5,.5,6.6],[53,2.2,9],[55,2,6],[49,1,8]]
    best=None
    for s in starts:
        res=minimize(fun,s,method="SLSQP",
                     bounds=[(19,55),(.5,3.5),(3,9)],
                     constraints=cons,
                     options={"ftol":1e-12,"maxiter":1500})
        if not res.success:
            continue
        p=predict(*res.x,recipe)
        feasible=(hull_margin(*res.x)>=-1e-6 and
                  max(p[r]-threshold for r in ["Coliforms","E. coli","Enterococci"])<=1e-6)
        if feasible and (best is None or res.fun<best.fun):
            best=res
    return best

def global_optimum(threshold,objective):
    candidates=[]
    for recipe in ["1S:3FW","2S:1FW","3S:1FW"]:
        res=optimise(recipe,threshold,objective)
        if res is not None:
            candidates.append((res.fun,recipe,res.x,predict(*res.x,recipe)))
    return min(candidates,key=lambda z:z[0])

reg=np.log10(1001)
jobs=[
    ("C1",reg,"min_temp"),
    ("C2",reg,"max_methane"),
    ("C3",2.5,"min_temp"),
    ("C4",2.5,"max_methane"),
]
rows=[]
for label,threshold,objective in jobs:
    _,recipe,z,p=global_optimum(threshold,objective)
    rows.append({
        "Scenario":label,
        "Temp_C":z[0],"OLR":z[1],"RT_days":z[2],"Recipe":recipe,
        "Coliforms":p["Coliforms"],"Ecoli":p["E. coli"],
        "Enterococci":p["Enterococci"],"Methane":p["Methane"],
        "Hull_margin":hull_margin(*z)
    })
out=pd.DataFrame(rows)
out_dir = ROOT / "outputs" / "tables"
out_dir.mkdir(parents=True, exist_ok=True)
out.to_csv(out_dir / "corrected_table3.csv", index=False)
print(out.to_string(index=False))
