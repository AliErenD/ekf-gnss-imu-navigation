function a = wrap_angle(a)
%WRAP_ANGLE Wrap angle(s) to [-pi, pi).
a = mod(a + pi, 2*pi) - pi;
end
