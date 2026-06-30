function plotSolutions(x,U,Uex)
plot(x,U,'r-d','LineWidth',1)
hold on
plot(x,Uex,'b-v','LineWidth',1)
legend('Solution EF','Solution exacte')
xlabel('x')
ylabel('u(x)')
grid on
end
