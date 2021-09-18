loadTS;
[sigma, Q, M] = FuncSignature(ts, 3, 1, 1);
x = 1:1:length(sigma);
stem(x,sigma);
xlabel('neuron label');
ylabel('population label');