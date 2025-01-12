% Load saved figures
big1=hgload('1big.fig');
small1=hgload('1small.fig');
big2=hgload('2big.fig');
small2=hgload('2small.fig');
big3=hgload('3big.fig');
small3=hgload('3small.fig');

% Prepare subplots
figure
h(1)=subplot(1,2,1);
h(2)=subplot(1,2,2);

figure
e(1)=subplot(1,2,1);
e(2)=subplot(1,2,2);

figure
f(1)=subplot(1,2,1);
f(2)=subplot(1,2,2);

% Paste figures on the subplots
copyobj(allchild(get(big1,'CurrentAxes')),h(1));
copyobj(allchild(get(small1,'CurrentAxes')),h(2));

copyobj(allchild(get(big2,'CurrentAxes')),e(1));
copyobj(allchild(get(small2,'CurrentAxes')),e(2));

copyobj(allchild(get(big3,'CurrentAxes')),f(1));
copyobj(allchild(get(small3,'CurrentAxes')),f(2));