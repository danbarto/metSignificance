import json

paraMC = [\
1.29,
1.193,
1.079,
1.135,
1.120,
-0.04,
0.6504
]


#paraMC = {\
#"a1":{'v':1.29,'eta':(0,0.8)},
#"a2":{'v':1.19,'eta':(0.8,1.3)},
#"a3":{'v':1.07,'eta':(1.3,1.9)},
#"a4":{'v':1.13,'eta':(1.9,2.5)},
#"a5":{'v':1.12,'eta':(2.5,100)},
#"N": {'v':-0.04, 'eta':(0,0)},
#"S": {'v':0.6504, 'eta':(0,0)}
#}

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

