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
	
###OBSOLETE, new versions have C's as lists corresponding to meas.
#IncludeSubgroup:= function(sub)
#local k;
#k := Length(sub);
#for i in [1..k] do
#	for j in [i..k] do
#		if sub[i]*sub[j] <> sub[j]*sub[i] then return false;
#		fi;
#	od;
#od;
#return true;
#end;

#AllSubgroups:= function(g)
#local ret;
#ret := Combinations(g);
##Removes the empty set
#Remove(ret,1);
#for h in [1..Length(g)] do
#	if IncludeSubgroup(ret[h]) then ;
#	else Remove(ret,h);
#	while IncludeSubgroup(ret[h]) <> true do Remove(ret,h);
#	od;
#	fi;
#od;
#return ret;
#end;

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

CompatibleSet:= function(g)
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

#Not needed with the state no longer being passed
#CompatibleSet:= function(state, meas)
#local compset;
#compset := AllSubgroups(meas);
#return compset;
#end;


#Part 2###############################################################

#calculates the empirical model - stored as a list of lists
#listed lists are probabilities ordered from (+++) -> (---)
#lists correspond to compatible sets of measurements

#probability uses projector with state vectors
#projectors are eigenvectors of each C-element

#For any (C1,...Cn), multiply each by each other if + or by 1- if -.

#EachSet := function(comps, meas)
#only 1 measurement is used
#if Length(comps) = 1 then return meas[comps[1]];
#local probs;

#Corresponds to res in original function.
Res := function(compset,meas,jp)
local res;
if jp[Length(jp)] = 0 then res := meas[compset[Length(jp)]];
else res := (IdentityMat(Length(meas[1])) - meas[compset[Length(jp)]]);
fi;
if Length(jp) > 1 then
	return Res(compset,meas,jp{[1..(Length(jp)-1)]}) * res;
else
	return res;
fi;
end;

#Separating this out.
BecauseLocal := function(comset,meas,jp)
local k;
k := Tuples([0,1],Length(comset));
return Res(comset,meas,k[jp]);
end;

Col := function(state, compset, meas)
local col;
col := [];
for j in [1..Length(Tuples([0,1],Length(compset)))] do
	#only multiply by the state at the end?
	#Print(BecauseLocal(compset,meas,j));
	Add(col,Trace(BecauseLocal(compset,meas,j)*state));
od;
return col;
end;

#Corresponds to col in original function.
#Col := function(state, compsets, meas, ip)
#local col;
#col := [];
#for j in [1..Length(Tuples([0,1],Length(compsets[ip])))] do
#	local temp;
#	temp := Res(compsets,meas,ip,k[j]);
#	Add(col,state*temp*state);
#od;
#return col;
#end;

###Further modded to separate res.
##Corresponds to col in original function.
#Col := function(state, compsets, meas, ip)
#local col;
#for j in [1..Tuples([0,1],Length(compsets[ip]))] do
#	local res;
#	res := ID();
#	for k in j do
#		local res1;
#		res1 := res;
#		if k = 0 then res := res1 * meas[compsets[ip]];
#		else res := res1 * (ID() - meas[compsets[ip]]);
#		fi;
#	od;
#	##Mult etc
#	Add(col,state*res*state);
#od;
#return col;
#end;

#Creates the empirical model
EmpMod := function(state, compsets, meas)
local grid;
grid := [];
#For each compatible set, performs this
for i in [1..Length(compsets)] do
	Add(grid,Col(state,compsets[i],meas));
od;
return grid;
end;

##Creates the empirical model
#EmpMod := function(state, compsets, meas)
#local grid;
##For each compatible set, performs this
#for i in [1..Length(compsets)] do
#	#res becomes the probability of the corresponding (?,...?).
#	local col;
#	for j in [1..Tuples([0,1],Length(compsets[i]))] do
#		local res;
#		res := ID();
#		for k in j do
#			local res1;
#			res1 := res;
#			if k = 0 then
#			res := res1 * meas[compsets[i]];
#			else
#			res := res1 * (ID() - meas[compsets[i]]);
#			fi;
#		od;
#		##Do the mult by the state/density matrix
#		Add(col,state*res*state);
#	od;
#	Add(grid,col);
#od;
#return grid;
#end;
	

#Part 3###############################################################

ConfirmZero := function(comp,empmodline,tuple)
local spot;
spot := Position(Tuples([0,1],Length(comp)),tuple{comp});
return empmodline[spot];
end;

#for every combination of [1,3] +/-'s, does every C have the EM undersaid column 0?
#comp is compatiblesets, empmod is the empirical model
ConCheck := function(meas,comp,empmod)
local flag;
for i in [1..Length(Tuples([0,1],Length(meas)))] do
	flag := false;
	for c in [1..Length(comp)] do
		if ConfirmZero(comp[c],empmod[c],Tuples([0,1],Length(meas))[i]) = 0 then
			flag := true;
			break;
		fi;
	od;
	if flag = false then
		Print("This is not contextual.\n");
		return;
		#return flag;
	fi;
od;
Print("This is contextual.\n");
#return true;
end;

#Overall##############################################################
#Setup for testing stage
nrow := function()
return [1/2,0,0,0,0,0,0,1/2];
end;

zrow := function()
return [0,0,0,0,0,0,0,0];
end;

pstate := function()
#return [nrow(),zrow(),zrow(),zrow(),zrow(),zrow(),zrow(),nrow()];
return [[1,0],[0,0]];
end;

zproj := function()
return [[1,1],[1,1]]*(1/2);
end;

pproj := function()
return [[1,-E(4)],[E(4),1]]*(1/2);
end;

xone := function()
return KroneckerProduct(zproj(),KroneckerProduct(ID(),ID()));
end;

xtwo := function()
return KroneckerProduct(ID(),KroneckerProduct(zproj(),ID()));
end;

xthr := function()
return KroneckerProduct(ID(),KroneckerProduct(ID(),zproj()));
end;

yone := function()
return KroneckerProduct(pproj(),KroneckerProduct(ID(),ID()));
end;

ytwo := function()
return KroneckerProduct(ID(),KroneckerProduct(pproj(),ID()));
end;

ythr := function()
return KroneckerProduct(ID(),KroneckerProduct(ID(),pproj()));
end;

pmeas := function()
return [xone(),xtwo(),xthr(),yone(),ytwo(),ythr()];
end;

automatic3 := function(comp,empmod)
local check;
ConCheck(pmeas(),comp,empmod);
#if check = true then Print("It is contextual.\n");
#else Print("It is not contextual.\n");
#fi;	
end;

automatic2 := function(comp)
local emp;
emp := EmpMod(pstate(),comp,pmeas());
automatic3(comp,emp);
end;

automatic := function()
local comps;
comps := CompatibleSet(pmeas());
automatic2(comps);
end;

performtest := function(state,measurements)
local com;
com := CompatibleSet(measurements);
ConCheck(measurements,com,EmpMod(state,com,measurements));
end;
