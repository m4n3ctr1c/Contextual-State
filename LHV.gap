#Utility##############################################################
#Store these as functions instead of global variables.
#Part 1###############################################################

#Tests for commutativity between all members of a given subgroup.
#sub is the subgroup that is being tested for inclusion
#meas is the set of measurements
#Returns false if any elements are not commutative, true if they are all comm.
IncludeSubgroup:= function(sub,meas)
local k;
k := Length(sub);
for i in [1..k] do
	for j in [i..k] do
		if meas[sub[i]]*meas[sub[j]] <> meas[sub[j]]*meas[sub[i]] then return false;
		fi;
	od;
od;
return true;
end;

#finds all subgroups of the measurement set and returns the commutative ones.
#g is the set of measurements
CompatibleSets:= function(g)
local ret;
ret := Combinations([1..Length(g)],LogInt(Length(g[1]),2));
for h in [1..Length(ret)] do
	if h <= Length(ret) then
	if IncludeSubgroup(ret[h],g) then ;
	else Remove(ret,h);
	while ((Length(ret) >= h) and (IncludeSubgroup(ret[h],g) <> true)) do Remove(ret,h);
	od;
	fi;
	else ;
	fi;
od;
return ret;
end;

#Part 2###############################################################

#Multiplies the measurements contained in a given compatible set.
#compset is the set of indices of compatible measurements
#meas is the set of measurements
#jp is the desired result of +/-'s
#Returns the product.
Res := function(compset,meas,jp)
local res;
#If the result of the measurement should be -, the measurement subtracted from the
#identity matrix instead.
if jp[Length(jp)] = 0 then res := meas[compset[Length(jp)]];
else res := (IdentityMat(Length(meas[1])) - meas[compset[Length(jp)]]);
fi;
#Multiplies from the back, since order does not matter.
if Length(jp) > 1 then
	return Res(compset,meas,jp{[1..(Length(jp)-1)]}) * res;
else
	return res;
fi;
end;

#Just a connecter due to limits on local variables.
BecauseLocal := function(comset,meas,jp)
local k;
k := Tuples([0,1],Length(comset));
return Res(comset,meas,k[jp]);
end;

Col := function(state, compset, meas)
local col;
col := [];
for j in [1..Length(Tuples([0,1],Length(compset)))] do
	Add(col,Trace(BecauseLocal(compset,meas,j)*state));
od;
return col;
end;

#Creates and returns the empirical model.
#state is the density matrix being tested
#compsets is the set of compatible groups
#meas is the set of measurements
EmpMod := function(state, compsets, meas)
local grid;
grid := [];
for i in [1..Length(compsets)] do
	Add(grid,Col(state,compsets[i],meas));
od;
return grid;
end;
	
#Part 3###############################################################

#GetRow := function(emprow,comp,tuples)
#local row;
#row := [];
#for j in [1..Length(tuples)] do
#	if emprow[Position(Tuples([0,1],Length(comp)),tuples[j]{comp})] = 0 then
#		Add(row,0);
#	else Add(row,1);
#	fi;
#od;
#return row;
#end;
GetRow := function(comp,j,tuples)
local row;
row := [];
for k in [1..Length(tuples)] do
	if tuples[k]{comp} = Tuples([0,1],Length(comp))[j] then Add(row,1);
	else Add(row,0);
	fi;
od;
return row;
end;

#BS := function(inc)
#local sum;
#sum := 0;
#for inci in [1..Length(inc)] do
#	sum := sum + inc[inci];
#od;
#return sum;
#end;

SoluteStuff := function(inc,svec)
Print(SolutionMat(inc,svec));
Print("\n");
Print(RankMat(inc));
Print("\n");
end;

CompareT := function(empmod,incidences)
local s;
s := [];
for j in [1..Length(empmod)] do
	for k in [1..Length(empmod[j])] do
		Add(s,empmod[j][k]);
	od;
od;
#Print(incidences);
#Print("\n");
#Print(s);
#Print("\n");
#Print(Length(s));
#Print("\n");
#Print(Length(incidences));
#Print(", ");
#Print(Length(incidences[1]));
#Print("\n");
#Print(BS(incidences[Length(incidences)]));
#Print("\n");
#Print(SolutionMat(TransposedMat(incidences),s));
#Print("\n");
SoluteStuff(TransposedMat(incidences),s);
end;

LHVTest := function(empmod,comps,meas)
local incidence;
incidence := [];
for i in [1..Length(empmod)] do
	for j in [1..Length(empmod[i])] do
	#Add(incidence,GetRow(empmod[i],comps[i],Tuples([0,1],meas)));
	Add(incidence,GetRow(comps[i],j,Tuples([0,1],meas)));
	od;
od;
##CREATING S AND SOLVING FOR T
#Print(empmod);
#Print("\n");
#Print(incidence);
#Print("\n");
CompareT(empmod,incidence);
end;

#Final Test####################################################################

performtest := function(state,measurements)
local com;
com := CompatibleSets(measurements);
LHVTest(EmpMod(state,com,measurements),com,Length(measurements));
end;

go := function(st,me)
local a;
a := Runtime();
performtest(st,me);
a := Runtime() - a;
Print(StringTime(a));
Print("\n");
end;
