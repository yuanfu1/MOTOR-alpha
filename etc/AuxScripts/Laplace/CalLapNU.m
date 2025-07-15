close all;
clear all;

% Parameters
h1 = 0;
h2 = 1;
h3 = 2;
h4 = 3;

a1 = (-2*(h2 + h3 + h4)/((h1 - h2)*(h1 - h3)*(h1 - h4)))
a2 = (2*(h1 + h3 + h4)/((h1 - h2)*(h2 - h3)*(h2 - h4)))
a3 = (2*(h1 + h2 + h4)/((h1 - h3)*(h3 - h2)*(h3 - h4)))
a4 = (2*(h1 + h2 + h3)/((h1 - h4)*(h4 - h2)*(h4 - h3)))

% 1st order derivative
h1 = -1;
h2 = 0;
h3 = 1;

a1 = -(h2+h3)/(h1-h2)/(h1-h3)
a2 = (h1+h3)/(h1-h2)/(h2-h3)
a3 = (h1+h2)/(h1-h3)/(h3-h2)

h1 = 0;
h2 = 1;
h3 = 2;

a1 = -(h2+h3)/(h1-h2)/(h1-h3)
a2 = (h1+h3)/(h1-h2)/(h2-h3)
a3 = (h1+h2)/(h1-h3)/(h3-h2)

h1 = -2;
h2 = -1;
h3 = 0;

a1 = -(h2+h3)/(h1-h2)/(h1-h3)
a2 = (h1+h3)/(h1-h2)/(h2-h3)
a3 = (h1+h2)/(h1-h3)/(h3-h2)




