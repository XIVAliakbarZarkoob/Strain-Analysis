

R = sort(R); hold on
plot([R(1) R(end)],[0.5 0.5],"r","LineWidth",3)
for DD = 0:150:5000
    
    LL = exp(-R.^2/DD^2);
    plot(R,LL,".","MarkerSize",14)
    title(sprintf("D = %g",DD))
    grid on
    pause

    if sum(LL > 0.5) > size(R,1)/2
        break
    end

end
title("Finding optimum value for D")
xlabel("R_{i}")
ylabel("L_{i}")
xlim([R(1) R(end)])



























