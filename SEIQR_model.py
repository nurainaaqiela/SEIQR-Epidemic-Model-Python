#!/usr/bin/env python
# coding: utf-8

# In[107]:


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# In[108]:


Lambda = 0.0002923
beta = 0.00003
mu = 0.0001077
theta = 0.00004
sigma = 0.1042308906
p = 0.5
gammaA = 0.0001
gammaS = 0.0001
gammaQ = 0.1
delta1 = 0.05
delta2 = 0.1


# In[109]:


def seiqr_model(t, y, beta):
    S, E, IA, IS, Q, R = y
    
    N = S + E + IA + IS + Q + R
    infection = beta * S * (IS + theta * IA) / N
    
    dSdt = Lambda - infection - mu*S
    dEdt = infection - (sigma + mu + delta1)*E
    dIAdt = p*sigma*E - (gammaA + mu)*IA
    dISdt = (1 - p)*sigma*E - (gammaS + mu + delta2)*IS
    dQdt = delta1*E + delta2*IS - (gammaQ + mu)*Q
    dRdt = gammaA*IA + gammaS*IS + gammaQ*Q - mu*R
    
    return [dSdt, dEdt, dIAdt, dISdt, dQdt, dRdt]


# In[110]:


N_equilibrium = Lambda / mu

S0 = N_equilibrium - 0.001
E0 = 0.0004
IA0 = 0.0003
IS0 = 0.0003
Q0 = 0
R0 = 0

initial_conditions = [S0, E0, IA0, IS0, Q0, R0]


# In[111]:


T = 180
t_span = (0, T)
t_eval = np.linspace(0, T, 1000)


# In[112]:


def compute_R0(beta_value):
    numerator = beta_value * Lambda * sigma
    denominator = mu * (sigma + mu + delta1)
    bracket = (p * theta) / (gammaA + mu) + (1 - p) / (gammaS + mu + delta2)
    return (numerator / denominator) * bracket


# In[113]:


# Baseline (R0 < 1)
beta_low = beta

# Increased transmission (force R0 > 1)
beta_high = 0.1


# In[114]:


R0_low = compute_R0(beta_low)
R0_high = compute_R0(beta_high)

print("R0 (Baseline) =", format(R0_low, ".6f"))
print("R0 (High Transmission) =", format(R0_high, ".6f"))


# In[115]:


# Compute proportionality constant
C = compute_R0(1)
print("Proportionality constant C =", C)

# Compute threshold beta where R0 = 1
beta_threshold = 1 / C
print("Beta threshold for R0 = 1:", beta_threshold)


# In[116]:


beta_high = beta_threshold * 5   

print("New beta_high =", beta_high)
print("New R0 =", compute_R0(beta_high))


# In[117]:


sol_low = solve_ivp(
    seiqr_model,
    t_span,
    initial_conditions,
    t_eval=t_eval,
    args=(beta_low,)  
)

sol_high = solve_ivp(
    seiqr_model,
    t_span,
    initial_conditions,
    t_eval=t_eval,
    args=(beta_high,)
)

# --- Extract symptomatic infected ---
IS_low = sol_low.y[3]
IS_high = sol_high.y[3]

# --- Peak detection ---
peak_index = np.argmax(IS_high)
peak_time = sol_high.t[peak_index]
peak_value = IS_high[peak_index]

plt.figure(figsize=(10,6))

plt.plot(sol_low.t, IS_low,
         label=f'R₀ < 1  (β = {beta_low:.5f})',
         linewidth=2)

plt.plot(sol_high.t, IS_high,
         label=f'R₀ > 1  (β = {beta_high:.3f})',
         linewidth=2)

# Peak marker 
plt.scatter(peak_time, peak_value,
            facecolor='orange',
            edgecolor='black',
            linewidth=1.5,
            s=120,
            zorder=5)

plt.annotate(
    f'Peak = {peak_value:.4f}\nTime = {peak_time:.1f} days',
    xy=(peak_time, peak_value),
    xytext=(peak_time - 40, peak_value + 0.003),   # moved higher
    fontsize=10,
    bbox=dict(boxstyle='round',
              facecolor='white',
              edgecolor='black',
              alpha=0.95),
    arrowprops=dict(arrowstyle='->')
)

plt.ylim(0, peak_value * 1.15)

plt.xlabel('Time (days)', fontsize=12)
plt.ylabel('Symptomatic Infected (IS)', fontsize=12)
plt.title('Threshold Behaviour of SEIQR Model', fontsize=16)

plt.legend(loc='upper left', fontsize=9)
plt.grid(alpha=0.2)
plt.tight_layout()
plt.show()


# In[118]:


def seiqr_model(t, y, beta, delta1_val, delta2_val):
    S, E, IA, IS, Q, R = y
    
    N = S + E + IA + IS + Q + R
    infection = beta * S * (IS + theta * IA) / N
    
    dSdt = Lambda - infection - mu*S
    dEdt = infection - (sigma + mu + delta1_val)*E
    dIAdt = p*sigma*E - (gammaA + mu)*IA
    dISdt = (1 - p)*sigma*E - (gammaS + mu + delta2_val)*IS
    dQdt = delta1_val*E + delta2_val*IS - (gammaQ + mu)*Q
    dRdt = gammaA*IA + gammaS*IS + gammaQ*Q - mu*R
    
    return [dSdt, dEdt, dIAdt, dISdt, dQdt, dRdt]


# In[119]:


beta_test = beta_high


# In[120]:


# With quarantine (original thesis values)
sol_with_Q = solve_ivp(
    seiqr_model,
    t_span,
    initial_conditions,
    t_eval=t_eval,
    args=(beta_test, delta1, delta2)
)

# Without quarantine (set isolation to zero)
sol_no_Q = solve_ivp(
    seiqr_model,
    t_span,
    initial_conditions,
    t_eval=t_eval,
    args=(beta_test, 0, 0)
)


# In[121]:


IS_with_Q = sol_with_Q.y[3]
IS_no_Q = sol_no_Q.y[3]


# In[122]:


# Peak detection
peak_with_Q = np.max(IS_with_Q)
peak_no_Q = np.max(IS_no_Q)

# Time of peaks
time_with_Q = sol_with_Q.t[np.argmax(IS_with_Q)]
time_no_Q = sol_no_Q.t[np.argmax(IS_no_Q)]

print("Peak WITH quarantine:", round(peak_with_Q, 4))
print("Peak WITHOUT quarantine:", round(peak_no_Q, 4))


# In[123]:


reduction_percentage = (1 - peak_with_Q / peak_no_Q) * 100

print("Percentage reduction due to quarantine:",
      round(reduction_percentage, 2), "%")


# In[106]:


# --- Peak detection ---
peak_with_Q = np.max(IS_with_Q)
peak_no_Q = np.max(IS_no_Q)

time_with_Q = sol_with_Q.t[np.argmax(IS_with_Q)]
time_no_Q = sol_no_Q.t[np.argmax(IS_no_Q)]

reduction_percentage = (1 - peak_with_Q / peak_no_Q) * 100

plt.figure(figsize=(10,6))

plt.plot(sol_with_Q.t, IS_with_Q,
         label='With Quarantine',
         linewidth=2)

plt.plot(sol_no_Q.t, IS_no_Q,
         label='Without Quarantine',
         linewidth=2)

# Mark peaks
plt.scatter(time_with_Q, peak_with_Q,
            color='blue', edgecolor='black',
            s=90, zorder=5)

plt.scatter(time_no_Q, peak_no_Q,
            color='orange', edgecolor='black',
            s=90, zorder=5)

# Annotation placement
plt.annotate(f'Peak (No Q)\n{peak_no_Q:.3f}',
             xy=(time_no_Q, peak_no_Q),
             xytext=(time_no_Q-30, peak_no_Q*0.9),
             arrowprops=dict(arrowstyle='->'),
             fontsize=9)

plt.annotate(f'Peak (With Q)\n{peak_with_Q:.3f}',
             xy=(time_with_Q, peak_with_Q),
             xytext=(time_with_Q-30, peak_with_Q*1.2),
             arrowprops=dict(arrowstyle='->'),
             fontsize=9)

# Move reduction box slightly left and centered vertically
plt.text(20, peak_no_Q*0.65,
         f'Peak Reduction = {reduction_percentage:.1f}%',
         fontsize=12,
         bbox=dict(facecolor='white', alpha=0.85))

plt.xlabel('Time (days)', fontsize=12)
plt.ylabel('Symptomatic Infected (IS)', fontsize=12)
plt.title('Impact of Quarantine on Disease Transmission', fontsize=16)

plt.legend(fontsize=10)
plt.grid(alpha=0.2)
plt.tight_layout()
plt.show()


# In[ ]:





# In[ ]:




