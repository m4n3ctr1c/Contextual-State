#Utility##############################################################
#Store these as functions instead of global variables.
XM:= function()
return [[0,1],[1,0]];
end;

YM:= function()
return [[0,E(4)],[-E(4),0]];
end;

ZM:= function()
return [[1,0],[0,-1]];
end;

ID:= function()
return [[1,0],[0,1]];
end;

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
ret := Combinations([1..Length(g)]);
#Removes the empty set
Remove(ret,1);
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

ConfirmZero := function(comp,empmodline,tuple)
local spot;
spot := Position(Tuples([0,1],Length(comp)),tuple{comp});
return empmodline[spot];
end;

PrintPlus := function(set)
Print("(");
if set[1] = 0 then Print("+");
else Print("-");
fi;
for l in [2..Length(set)] do
	Print(",");
	if set[l] = 0 then Print("+");
	else Print("-");
	fi;
od;
Print("): ");
end;

#for every combination of [1,3] +/-'s, does every C have the EM undersaid column 0?
#comp is compatiblesets, empmod is the empirical model
ConCheck := function(meas,comp,empmod)
local flag;
flag := true;
for i in [1..Length(Tuples([0,1],meas))] do
	PrintPlus(Tuples([0,1],meas)[i]);
	for c in [1..Length(comp)] do
		if ConfirmZero(comp[c],empmod[c],Tuples([0,1],meas)[i]) = 0 then
			Print(comp[c]);
			Print("\n");
			break;
		fi;
		if c = Length(comp) then
			flag := false;
			Print("None\n");
		fi;
	od;
od;
if flag = true then Print("This is contextual.\n");
else Print("This is not contextual.\n");
fi;
end;

#Final Test#########################################################################

performtest := function(state,measurements)
local com;
com := CompatibleSets(measurements);
ConCheck(Length(measurements),com,EmpMod(state,com,measurements));
end;
