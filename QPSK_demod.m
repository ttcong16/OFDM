function bits = QPSK_demod(rx)
N = length(rx);
bits = zeros(2*N,1);
for k = 1:N
    re = real(rx(k)); im = imag(rx(k));
    if re>=0 && im>=0, bits(2*k-1:2*k)=[0;0];
    elseif re<0 && im>=0, bits(2*k-1:2*k)=[0;1];
    elseif re<0 && im<0, bits(2*k-1:2*k)=[1;1];
    else bits(2*k-1:2*k)=[1;0];
    end
end
end
