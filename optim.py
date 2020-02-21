# find optimal fuzzy solution for three equations in two unknowns
# using sum of squares

bignum = 1000000000000

# coefficients by call

coeff = {}
coeff["11"] = coeff["22"] = coeff["33"] = 2.0
coeff["12"] = coeff["24"] = coeff["36"] = 1.0
coeff["13"] = 2.0/3.0
coeff["23"] = 4.0/3.0
coeff["14"] = 1.0/2.0
coeff["34"] = 6.0/4.0
coeff["15"] = 2.0/5.0
coeff["25"] = 4.0/5.0
coeff["35"] = 6.0/5.0
coeff["16"] = 1.0/3.0

class eqn:
  def __init__(self,eqnstring):
    # parse eqnstring into variables and call
    self.parameters = []
    e = eqnstring.split("_")
    self.call = "12"
    for item in e:
      if item == "x1":
        self.parameters.append(0)
      elif item == "x2":
        self.parameters.append(1)
      elif item == "x3":
        self.parameters.append(2)
      elif item == "x4":
        self.parameters.append(3)
      else:    # this term is optional
        self.call = item

    assert len(self.parameters) > 0 and len(self.parameters) <= 4

    # figure out the right coefficient
    self.coeff = coeff[self.call]

  def __call__(self,inputs):
    inputlist = []
    for index in self.parameters:
      inputlist.append(inputs[index])
    result = sum(inputlist) * self.coeff
    return result

class datum:
  def __init__(self):
    self.weights = []
    self.results = []
    self.equations = []
    self.eqtypes = []
  def add(self,result,weight,eqtype):
    self.results.append(result)
    self.weights.append(weight)
    self.eqtypes.append(eqtype)
    self.equations.append(eqn(eqtype))
    
def which_variables_present(equationlist):
  x1_mentioned = False
  x2_mentioned = False
  x3_mentioned = False
  x4_mentioned = False
  for item in equationlist:
    if "x1" in item:  x1_mentioned = True
    if "x2" in item:  x2_mentioned = True
    if "x3" in item:  x3_mentioned = True
    if "x4" in item:  x4_mentioned = True
  mentions = []
  if x1_mentioned:  mentions.append("x1")
  if x2_mentioned:  mentions.append("x2")
  if x3_mentioned:  mentions.append("x3")
  if x4_mentioned:  mentions.append("x4")
  return mentions
  
def getscore(x1,x2,x3,x4,equations,weights,results):
  score = 0.0
  inputlist = [x1,x2,x3,x4]
  for eq,weight,result in zip(equations,weights,results):
    score += (eq(inputlist) - result)**2.0 * weight
  return score

def getprefix(mentionlist):
  prefixes = []
  if "x1" in mentionlist:
    prefixes.append("1")
  if "x2" in mentionlist:
    prefixes.append("2")
  if "x3" in mentionlist:
    prefixes.append("3")
  if "x4" in mentionlist:
    prefixes.append("4")
  return prefixes


###########################################################################
# program

import sys
if len(sys.argv) != 2:
  print "Usage:  python optim2.py datafile"
  exit()

datadict = {}
datadict_by_call = {}
for line in open(sys.argv[1],"r"):
  if line == "\n":  continue
  if line.startswith("#"):  continue     # commented out lines
  line = line.split()
  if len(line) == 4:
    code,result,weight,eqtype = line
    call = "1,1"
  else:
    code,result,weight,eqtype,call = line
  result = float(result)
  weight = float(weight)
  if code not in datadict:
    datadict[code] = datum()
  if code not in datadict_by_call:
    datadict_by_call[code] = {}
  if call not in datadict_by_call[code]:
    datadict_by_call[code][call] = datum()
  datadict[code].add(result,weight,eqtype)
  datadict_by_call[code][call].add(result,weight,eqtype)

cf_range = range(1,101)

keylist = datadict.keys()
keylist.sort()

print
print "OVERALL"
print "-"*40
for code in keylist:
  mydata = datadict[code]
  results = mydata.results
  weights = mydata.weights
  equations = mydata.equations
  eqtypes = mydata.eqtypes

  mentions = which_variables_present(eqtypes)
  minscore = bignum
  minx1 = None
  minx2 = None
  minx3 = None
  minx4 = None

  # One dimensional search
  if len(mentions) == 1:
    for x in cf_range:
      score = getscore(x,0,0,0,equations,weights,results)
      if score < minscore:
        minx = x
        minscore = score
    if minscore == bignum: 
      print "No maximum found for ",code
    outstring = code + ": " + code + getprefix(mentions)[0] + "=" + str(minx)
    outstring += " score " + str(minscore)
    print outstring
    continue

  # Two dimensional search
  if len(mentions) == 2:
    x3 = 0
    x4 = 0
    for x1 in cf_range:
      for x2 in cf_range:
        if x1 + x2 > 100:  continue
        score = getscore(x1,x2,x3,x4,equations,weights,results)
        if score < minscore:
          minx1 = x1
          minx2 = x2
          minscore = score
    if minscore == bignum: 
      print "No maximum found for ",code
    outstring = code + ": " + code + getprefix(mentions)[0] + "=" + str(minx1) 
    outstring += ", " + code + getprefix(mentions)[1] + "=" + str(minx2) 
    outstring += " score " + str(minscore)
    print outstring
    continue

  # Three dimensional search
  if len(mentions) == 3:
    x4 = 0
    for x1 in cf_range:
      for x2 in cf_range:
        for x3 in cf_range:
          if x1 + x2 + x3 > 100:  continue
          score = getscore(x1,x2,x3,x4,equations,weights,results)
          if score < minscore:
            minx1 = x1
            minx2 = x2
            minx3 = x3
            minscore = score
    if minscore == bignum: 
      print "No maximum found for ",code
    outstring = code + ": " + code + getprefix(mentions)[0] + "=" + str(minx1) 
    outstring += ", " + code + getprefix(mentions)[1] + "=" + str(minx2) 
    outstring += ", " + code + getprefix(mentions)[2] + "="  + str(minx3)
    outstring += " score " + str(minscore)
    print outstring

  # Four dimensional search!
  if len(mentions) == 4:
    for x1 in cf_range:
      for x2 in cf_range:
        for x3 in cf_range:
          for x4 in cf_range:
            if x1 + x2 + x3 + x4 > 100:  continue
            score = getscore(x1,x2,x3,x4,equations,weights,results)
            if score < minscore:
              minx1 = x1
              minx2 = x2
              minx3 = x3
              minx4 = x4
              minscore = score
    if minscore == bignum:
      print "No maximum found for ",code
    outstring = code + ": " + code + getprefix(mentions)[0] + "=" + str(minx1) 
    outstring += ", " + code + getprefix(mentions)[1] + "=" + str(minx2) 
    outstring += ", " + code + getprefix(mentions)[2] + "="  + str(minx3)
    outstring += ", " + code + getprefix(mentions)[3] + "="  + str(minx4)
    outstring += " score " + str(minscore)
    print outstring

# same analysis broken out by call
print
print "BY CALL"
print "-"*40
outputs = {}
for code in keylist:
  if code not in outputs:
    outputs[code] = []
  for call in datadict_by_call[code]:
    mydata = datadict_by_call[code][call]
    results = mydata.results
    weights = mydata.weights
    equations = mydata.equations
    eqtypes = mydata.eqtypes

    mentions = which_variables_present(eqtypes)
    minscore = bignum
    minx1 = None
    minx2 = None
    minx3 = None
    minx4 = None

    # One dimensional search
    if len(mentions) == 1:
      for x in cf_range:
        score = getscore(x,0,0,0,equations,weights,results)
        if score < minscore:
          minx = x
          minscore = score
      if minscore == bignum: 
        print "No maximum found for ",code
      outstring = code + ": " + code + getprefix(mentions)[0] + "=" + str(minx) + " " + call
      outputs[code].append(outstring)
      continue

    # Two dimensional search
    if len(mentions) == 2:
      x3 = 0
      x4 = 0
      for x1 in cf_range:
        for x2 in cf_range:
          if x1 + x2 > 100:  continue
          score = getscore(x1,x2,x3,x4,equations,weights,results)
          if score < minscore:
            minx1 = x1
            minx2 = x2
            minscore = score
      if minscore == bignum: 
        print "No maximum found for ",code
      outstring = code + ": " + code + getprefix(mentions)[0] + "=" + str(minx1)
      outstring += ", " + code + getprefix(mentions)[1] + "=" + str(minx2) + " " + call 
      outputs[code].append(outstring)
      continue

    # Three dimensional search
    if len(mentions) == 3:
      x4 = 0
      for x1 in cf_range:
        for x2 in cf_range:
          for x3 in cf_range:
            if x1 + x2 + x3 > 100:  continue
            score = getscore(x1,x2,x3,x4,equations,weights,results)
            if score < minscore:
              minx1 = x1
              minx2 = x2
              minx3 = x3
              minscore = score
      if minscore == bignum: 
        print "No maximum found for ",code
      outstring = code + ": " + code + getprefix(mentions)[0] + "=" + str(minx1) 
      outstring += ", " + code + getprefix(mentions)[1] + "=" + str(minx2) 
      outstring += ", " + code + getprefix(mentions)[2] + "="  + str(minx3) + " " + call
      outputs[code].append(outstring)

    # Four dimensional search!
    if len(mentions) == 4:
      for x1 in cf_range:
        for x2 in cf_range:
          for x3 in cf_range:
            for x4 in cf_range:
              if x1 + x2 + x3 + x4 > 100:  continue
              score = getscore(x1,x2,x3,x4,equations,weights,results)
              if score < minscore:
                minx1 = x1
                minx2 = x2
                minx3 = x3
                minx4 = x4
                minscore = score
      if minscore == bignum:
        print "No maximum found for ",code
      outstring = code + ": " + code + getprefix(mentions)[0] + "=" + str(minx1) 
      outstring += ", " + code + getprefix(mentions)[1] + "=" + str(minx2) 
      outstring += ", " + code + getprefix(mentions)[2] + "="  + str(minx3)
      outstring += ", " + code + getprefix(mentions)[3] + "="  + str(minx4)
      outputs[code].append(outstring)


for code in keylist:
  for entry in outputs[code]:
    print entry
  print "-"*40
