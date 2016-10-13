import json

paraMC = {\
"a1":1.29,
"a2":1.19,
"a3":1.07,
"a4":1.13,
"a5":1.12,
"N":-0.04,
"S":0.6504
}

with open('../data/paraMC.txt', 'w') as outfile:
  json.dump(paraMC, outfile)

paraData = {\
"a1":1.26,
"a2":1.14,
"a3":1.13,
"a4":1.13,
"a5":1.06,
"N":-3.3,
"S":0.5961
}

with open('../data/paraData.txt', 'w') as outfile:
  json.dump(paraData, outfile)

