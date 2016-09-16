function out = fett(figNr)

% macht alles fetter (Linien, Text etc).

figure(figNr);

lw = 4;   %line width
fs = 18;  %fontsize   

h = findobj(gca,'Type','line');
set(h,'LineWidth',lw);

set(gca,'FontSize',fs);
set(gca,'LineWidth',lw);
set(get(gca,'title'),'Fontsize',fs);  
set(get(gca,'xlabel'),'Fontsize',fs); 
set(get(gca,'ylabel'),'Fontsize',fs);   
set(legend,'Fontsize',fs);   
