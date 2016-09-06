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

Poss := function(ga,compset,tuple,size)
local tuples;
tuples := Tuples([0,1],size);
for o in [1..Length(tuples)] do
	if tuples[o]{compset} = tuple then
		if ga[o] = 1 then
	#Print(o);
	#Print(": ");
			return true;
		fi;
	fi;
od;
return false;
end;

SetZeros := function(ga,compset,tuple,size)
local tuples;
tuples := Tuples([0,1],size);
for l in [1..Length(ga)] do
	if tuples[l]{compset} = tuple then
		ga[l] := 0;
	fi;
od;
return ga;
end;

HardyTest := function(empmod,comps,meas)
local ga;
ga := [];
for i in [1..2^meas] do
	Add(ga,1);
od;
for j in [1..Length(empmod)] do
	for k in [1..Length(empmod[j])] do
		if empmod[j][k] = 0 then
			ga := SetZeros(ga,comps[j],Tuples([0,1],Length(comps[j]))[k],meas);
		fi;
	od;
od;
#Print(empmod);
#Print("\n");
#Print(ga);
#Print("\n");
for m in [1..Length(empmod)] do
	for n in [1..Length(empmod[m])] do
		if empmod[m][n] > 0 then
			if Poss(ga,comps[m],Tuples([0,1],Length(comps[m]))[n],meas) = false then
				#false means that it IS hardy
				Print("It is Hardy.\n");
				return;
#else
	#Print(m);
	#Print(",");
	#Print(n);
	#Print("\n");
			fi;
		fi;
	od;
od;
Print("It is not Hardy.\n");
end;

#Final Test####################################################################

performtest := function(state,measurements)
local com;
com := CompatibleSets(measurements);
HardyTest(EmpMod(state,com,measurements),com,Length(measurements));
end;

go := function(st,me)
local a;
a := Runtime();
performtest(st,me);
a := Runtime() - a;
Print(StringTime(a));
Print("\n");
end;
