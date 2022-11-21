neuros = [96, 97, 98, 99, 100, 110, 120, 150, 200, 500];
trans_time = [5000, 360, 240, 210, 190, 120, 90, 50, 33, 15];
plot(neuros, trans_time, 'b.-')
hold on;
plot([96, 96], [0, 5000], 'r--')
hold off;
title("Time until transmission as a function of neurotransmitters released")
xlabel("Neurotransmitters released into system")
ylabel("Time [ns]")
ylim([0 500])