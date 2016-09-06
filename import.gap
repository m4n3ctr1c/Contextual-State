ID := function(j)
return IdentityMat(2^j);
end;

#kth x, s operations
xg := function(k,s)
local zp;
zp := [[1,1],[1,1]]/2;
if k = 1 then
	if s = 1 then return zp;
	else return KroneckerProduct(zp,ID(s-1));
	fi;
elif k = s then return KroneckerProduct(ID(s-1),zp);
else return KroneckerProduct(ID(k-1),KroneckerProduct(zp,ID(s-k)));
fi;
end;

#kth y, s operations
yg := function(k,s)
local pp;
pp := [[1,-E(4)],[E(4),1]]/2;
if k = 1 then
	if s = 1 then return pp;
	else return KroneckerProduct(pp,ID(s-1));
	fi;
elif k = s then return KroneckerProduct(ID(s-1),pp);
else return KroneckerProduct(ID(k-1),KroneckerProduct(pp,ID(s-k)));
fi;
end;

#kth z, s operations
zg := function(k,s)
local pp;
pp := [[1,0],[0,0]];
if k = 1 then
	if s = 1 then return pp;
	else return KroneckerProduct(pp,ID(s-1));
	fi;
elif k = s then return KroneckerProduct(ID(s-1),pp);
else return KroneckerProduct(ID(k-1),KroneckerProduct(pp,ID(s-k)));
fi;
end;

#kth h, s operations
hg := function(k,s)
local hp;
hp := [[2+Sqrt(2),Sqrt(2)],[Sqrt(2),2-Sqrt(2)]]/4;
if k = 1 then
	if s = 1 then return hp;
	else return KroneckerProduct(hp,ID(s-1));
	fi;
elif k = s then return KroneckerProduct(ID(s-1),hp);
else return KroneckerProduct(ID(k-1),KroneckerProduct(hp,ID(s-k)));
fi;
end;

gg := function(obs,s)
local rets;
rets := [KroneckerProduct(obs,ID(s-1))];
for i in [2..s-1] do
	Add(rets,KroneckerProduct(ID(i-1),KroneckerProduct(obs,ID(s-i))));
od;
Add(rets,KroneckerProduct(ID(s-1),obs));
return rets;
end;
