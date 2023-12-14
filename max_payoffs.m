function[u1, u2] = max_payoffs(R1, R2, yo, zo, yn, zn)

if zo > zn
    u1 = zo*(-R1(2,2) + R1(2,1)) + R1(2,2);
else
    u1 = zo * R1(1,1);
end

if yo > yn
    u2 = yo * (R2(1,1) - R2(2,1)) + R2(2,1);
else
    u2 = 0;
end