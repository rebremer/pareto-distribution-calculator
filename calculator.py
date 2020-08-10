# For more information on Pareto distribution and its properties, refer to
#
# https://en.wikipedia.org/wiki/Pareto_distribution
# https://wiki.socr.umich.edu/index.php/AP_Statistics_Curriculum_2007_Pareto
# https://math.stackexchange.com/questions/2909207/continuous-pareto-distribution-intuition
# http://www.eclecticon.info/index_htm_files/Probability%20&%20Statistics%20-%20Expectation%20Variance%20Skew%20Kurtosis.pdf
# https://arxiv.org/abs/2001.10488

# Parameters
alpha = 5
x_max = pow(10,50) # Pareto distribution is really fat tailed for low values of alpha (close to one), so a large value shall be taken to compute E[x] accurately

# Variables
x_m = 1000
x = x_m
p_sum = 0
e_x = 0
e_x2 = 0
e_x3 = 0
e_x4 = 0
var_x = 0
iterator = 1

# Compute Expected Value E[x] and moments E[x^2], E[x^3] and E[x^4]
while x < x_max:

    p_x = (alpha * pow(x_m,alpha))/(pow(x,alpha + 1))
    p_sum += p_x * iterator
    e_x += x * p_x * iterator
    e_x2 += pow(x,2) * p_x * iterator
    e_x3 += pow(x,3) * p_x * iterator
    e_x4 += pow(x,4) * p_x * iterator

    x += 1 * iterator
    
    # Step size shall be multiplied by 10 for every 10^n iteration
    if x % (x_m * iterator * 10) == 0:
        iterator *= 10
        #print("iteration: " + str(x) + ", p_x: " + str(p_x) + ", expected value: " + str(e_x))

# Compute Var[x] and Kur[x] using moments
var_x = e_x2 - pow(e_x,2) # between alpha 4 and 8 computation is pretty accurate
kurtosis_x = (e_x4 - 4 * e_x * e_x3 + 6 * pow(e_x,2) * e_x2 - 3 * pow(e_x,4))/(pow(var_x,2)) -3 # between alpha 4 and 8 computation is pretty accurate

# Calculate E[x], Var[x] and Kur[x] using definitions for standard Pareto definition
if alpha > 1:
    e_xcalc = (alpha * x_m)/(alpha-1)
else:
    e_xcalc = "n/a for alpha =< 1"

if alpha > 2:
    var_xcalc = (alpha * pow(x_m,2))/(pow((alpha-1),2)*(alpha - 2))
else:
    var_xcalc = "n/a for alpha =< 2"

if alpha > 4:
    kur_xcalc = 6*(pow(alpha,3) + pow(alpha,2) - 6 * alpha - 2)/ (alpha * (alpha - 3) * (alpha - 4))
else:
    kur_xcalc = "n/a for alpha =< 4"

print("sum p_x: " + str(p_sum) + ", e_x_computed: " + str(e_x) + ", e_x_calculus: " + str(e_xcalc))
print("var_x_computed: " + str(var_x) + ", var_x_calculus: " + str(var_xcalc)) 
print("kur_computed: " + str(kurtosis_x) + ", kur_xcalc: " + str(kur_xcalc))
