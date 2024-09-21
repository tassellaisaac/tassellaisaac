Define the transition matrix for the Markov Chain.
Simulate the progression of a single patient (or multiple patients) through the states of the Markov Chain.
Plot or analyze the resulting simulations.

1. Define the States:
We have five states in the system:

Healthy (H): 1
Sick (S0): 2
Hospital (SH): 3
Dead (D): 4
2. Define the Transition Matrix:
We can represent the transitions between these states using a matrix, where each entry (i, j) represents the probability of moving from state i to state j.

H -> S0 = 0.2
H -> SH = 0.1
H -> H = ? (This will be 1 - (0.2 + 0.1) = 0.7)
S0 -> H = 0.2
S0 -> SH = 0.1
S0 -> D = 0.01
S0 -> S0 = ? (This will be 1 - (0.2 + 0.1 + 0.01) = 0.69)
SH -> H = 0.1
SH -> D = 0.2
SH -> SH = ? (This will be 1 - (0.1 + 0.2) = 0.7)
D -> D = 1

# Set seed for reproducibility
set.seed(123)

# Define states
states <- c("Healthy", "Sick", "Hospital", "Dead")

# Define the transition matrix
transition_matrix <- matrix(c(
  0.7, 0.2, 0.1, 0,     # From Healthy
  0.2, 0.69, 0.1, 0.01, # From Sick
  0.1, 0, 0.7, 0.2,     # From Hospital
  0, 0, 0, 1            # From Dead
), nrow = 4, byrow = TRUE)

# Set row and column names for readability
rownames(transition_matrix) <- states
colnames(transition_matrix) <- states

# Function to simulate Markov Chain
simulate_markov_chain <- function(initial_state, steps) {
  current_state <- initial_state
  state_sequence <- c(current_state)
  
  for (i in 1:steps) {
    # Get probabilities for the current state
    current_state_probs <- transition_matrix[current_state, ]
    
    # Sample the next state based on the transition probabilities
    next_state <- sample(1:4, 1, prob = current_state_probs)
    
    # Append the new state to the sequence
    state_sequence <- c(state_sequence, next_state)
    
    # Update the current state
    current_state <- next_state
  }
  
  # Return the sequence of states
  return(state_sequence)
}

# Simulate for a single patient starting from Healthy (state 1)
simulated_chain <- simulate_markov_chain(1, 20)

# Convert numeric states to state names
state_names_sequence <- states[simulated_chain]

# Print the sequence of states
print(state_names_sequence)

You can simulate multiple patients or run longer simulations by increasing the number of steps. The current example simulates 20 steps for a single patient starting in the Healthy state.

# Plot the state transitions
plot(state_names_sequence, type = "b", main = "Markov Chain Simulation", 
     xlab = "Step", ylab = "State", col = "blue", pch = 19, xaxt = "n")

# Set custom x-axis labels
axis(1, at = 1:length(state_names_sequence), labels = 1:length(state_names_sequence))

You can change the number of steps and the initial state by modifying the inputs to simulate_markov_chain().
You can simulate multiple patients by calling the function multiple times and analyzing the results collectively.

Simplified Three-State Model:
State 1 (H): Healthy
State 2 (S0): Sick
State 3 (D): Dead

Healthy (H):

𝑃
(
𝐻
→
𝑆
0
)
=
0.2
P(H→S0)=0.2
𝑃
(
𝐻
→
𝐻
)
=
𝑎
P(H→H)=a (unknown)
𝑃
(
𝐻
→
𝐷
)
=
1
0
−
3
P(H→D)=10 
−3
 
The total probability must sum to 1, so 
𝑎 = 1 −(0.2 + 10 − 3) = 0.799
a = 1−(0.2+10 −3) =0.799
Sick (S0): 𝑃(𝑆0→𝐻) = 0.2
P(S0→H)=0.2𝑃(𝑆0→𝑆0)=𝑏
P(S0→S0)=b (unknown)
𝑃(𝑆0→𝐷)=0.01
P(S0→D)=0.01
The total probability must sum to 1, so 𝑏=1−(0.2+0.01)=0.79
b=1−(0.2+0.01)=0.79
Dead (D):𝑃(𝐷→𝐷)=1
P(D→D)=1 (no other transitions are possible from the dead state).

# Set seed for reproducibility so that each run gives the same result.
set.seed(42)

# Define the number of steps to simulate in the Markov Chain.
n_steps <- 10  # You can change this to any number for longer simulations.

# Initialize a vector to store the sequence of states during the simulation. 
# Start with the first state (Healthy), which is state 1.
state_sequence <- numeric(n_steps)  
state_sequence[1] <- 1  # Start from the Healthy state.

# Define the transition matrix for the Markov Chain.
# Rows correspond to the current state, and columns correspond to the next state.
transition_matrix <- matrix(c(
  # From Healthy (H): P(H -> H), P(H -> S0), P(H -> D)
  0.799, 0.2, 0.001, 
  # From Sick (S0): P(S0 -> H), P(S0 -> S0), P(S0 -> D)
  0.2, 0.79, 0.01, 
  # From Dead (D): P(D -> H), P(D -> S0), P(D -> D)
  0, 0, 1
), nrow = 3, byrow = TRUE)

# Loop through the number of steps to simulate the Markov Chain.
for (i in 2:n_steps) {
  # Get the current state from the state sequence.
  current_state <- state_sequence[i - 1]
  
  # Sample the next state based on the current state's transition probabilities.
  next_state <- sample(1:3, size = 1, prob = transition_matrix[current_state, ])
  
  # Store the next state in the state sequence.
  state_sequence[i] <- next_state
}

# Print the sequence of states as they evolve over time.
print(state_sequence)

# Convert numeric states to meaningful labels for readability.
states <- c("Healthy", "Sick", "Dead")
state_labels <- states[state_sequence]

# Print the sequence with state labels.
print(state_labels)

set.seed(42)
state_sequence <- numeric(n_steps)
state_sequence[1] <- 1  # Start from Healthy
transition_matrix <- matrix(c(
  0.799, 0.2, 0.001,   # Healthy: H -> H, H -> S0, H -> D
  0.2, 0.79, 0.01,     # Sick: S0 -> H, S0 -> S0, S0 -> D
  0, 0, 1              # Dead: D -> D
), nrow = 3, byrow = TRUE)
for (i in 2:n_steps) {
  current_state <- state_sequence[i - 1]
  next_state <- sample(1:3, size = 1, prob = transition_matrix[current_state, ])
  state_sequence[i] <- next_state
}
states <- c("Healthy", "Sick", "Dead")
state_labels <- states[state_sequence]
print(state_labels)

a=0.799: The probability of staying in the Healthy state.
𝑏=0.79
b=0.79: The probability of staying in the Sick (S0) state.

This R script simulates a three-state Markov Chain representing the progression of a person’s health from Healthy to Sick to Dead. The probabilities for each transition are defined in the transition matrix, and the simulation evolves over a specified number of steps.

State 1 (H): Healthy
State 2 (S0): Sick (managed at home)
State 3 (SH): Hospitalized
State 4 (D): Dead

set.seed(42)

# Define the transition probabilities for the 4-state Markov Chain
a <- 0.7  # Probability of staying Healthy (H)
b <- 0.69 # Probability of staying Sick (S0)

# Number of days to simulate
n_days <- 400

# Define the 4x4 transition matrix
# From state H, S0, SH, D
transition_matrix <- matrix(c(
  # H: Healthy
  a, 0.2, 0.1, 0.001,
  # S0: Sick at home
  0.2, b, 0.1, 0.01,
  # SH: Hospitalized
  0.1, 0, 0.7, 0.2,
  # D: Dead
  0, 0, 0, 1
), nrow = 4, ncol = 4, byrow = TRUE)

# Initialize the patient's starting state (Healthy)
state <- 1  # Starting at Healthy (H)

# Initialize a vector to keep track of the patient's state for each day
patient_record <- rep(0, n_days)

# Simulate the Markov Chain over n_days
for (day in 1:n_days) {
  # Get the transition probabilities for the current state
  pr <- transition_matrix[state, ]
  
  # Sample the next state based on the current state's transition probabilities
  state <- sample(c(1:4), size = 1, prob = pr)
  
  # Record the state for the current day
  patient_record[day] <- state
}

# Plot the resulting patient state record over time
plot(1:n_days, patient_record, type = "l", col = "blue", lwd = 2,
     main = "4-State Markov Chain: Patient Health Over Time",
     xlab = "Days", ylab = "State", yaxt = "n")

# Add labels to the y-axis
axis(2, at = 1:4, labels = c("Healthy (H)", "Sick (S0)", "Hospital (SH)", "Dead (D)"))

Extended to a 4x4 matrix to include the hospitalization state (SH).
Probabilities sum to 1 for each row. For example, in the Healthy (H) state, the patient can:
Stay healthy with probability a = 0.7
Become sick at home with probability 0.2
Be hospitalized with probability 0.1
Die with probability 0.001

The patient's state is simulated for 400 days using the new transition matrix.
At each step, the next state is chosen based on the current state's transition probabilities.

The y-axis represents the patient's state, where:
1 = Healthy (H)
2 = Sick (S0)
3 = Hospitalized (SH)
4 = Dead (D)
The resulting plot shows how the patient's health transitions over time.

Healthy Periods: The patient may remain healthy for several days (state 1) due to the high probability of staying healthy.
Sick and Hospitalization Cycles: After a period of being healthy, the patient may become sick (state 2) and eventually be hospitalized (state 3), depending on the transition probabilities.
Death State: Once the patient transitions to the Dead (D) state, they remain in that state, as indicated by a horizontal line at state 4.

The simulation will likely show the patient cycling between Healthy (H), Sick (S0), and Hospitalized (SH) for a while. However, at some point, the patient may transition to the Dead (D) state and stay there indefinitely (due to the absorbing nature of death in the model).
The patient may stay healthy for many days at a time but can eventually become sick or hospitalized, depending on the transition probabilities.

set.seed(42)

# Define the transition probabilities for the 4-state Markov Chain
a <- 0.7  # Probability of staying Healthy (H)
b <- 0.69 # Probability of staying Sick (S0)

# Number of days to simulate
n_days <- 400

# Define the 4x4 transition matrix
transition_matrix <- matrix(c(
  # H: Healthy
  a, 0.2, 0.1, 0.001,
  # S0: Sick at home
  0.2, b, 0.1, 0.01,
  # SH: Hospitalized
  0.1, 0, 0.7, 0.2,
  # D: Dead
  0, 0, 0, 1
), nrow = 4, ncol = 4, byrow = TRUE)

# Function to simulate the Markov Chain for a single patient
simulate_patient <- function(n_days, transition_matrix) {
  state <- 1  # Start at Healthy (state 1)
  patient_record <- rep(0, n_days)
  
  for (day in 1:n_days) {
    pr <- transition_matrix[state, ]
    state <- sample(1:4, size = 1, prob = pr)
    patient_record[day] <- state
  }
  
  return(patient_record)
}

# Simulate 1000 patient records
n_patients <- 1000
all_patient_records <- matrix(0, nrow = n_patients, ncol = n_days)

for (i in 1:n_patients) {
  all_patient_records[i, ] <- simulate_patient(n_days, transition_matrix)
}

# Calculate the number of days spent in each state for all patients
days_in_states <- matrix(0, nrow = n_patients, ncol = 4)

for (i in 1:n_patients) {
  days_in_states[i, ] <- table(factor(all_patient_records[i, ], levels = 1:4))
}

# Sum the total number of days spent in each state across all patients
total_days_in_states <- colSums(days_in_states)

# Plot the distribution of the number of days spent in each state
barplot(total_days_in_states, names.arg = c("Healthy (H)", "Sick (S0)", "Hospital (SH)", "Dead (D)"),
        col = c("green", "yellow", "orange", "red"), 
        main = "Distribution of Days Spent in Each State (1000 Patients)",
        ylab = "Total Number of Days", xlab = "State")

Healthy (H): Since the probability of staying healthy is relatively high (a = 0.7), we expect most patients to spend a significant number of days in the Healthy state.
Sick (S0): Patients transition to the Sick state with a moderate probability of 0.2, and can stay sick with probability b = 0.69. This state will likely show fewer days than the Healthy state but still a substantial amount.
Hospital (SH): The probability of transitioning to Hospital is lower compared to becoming sick, and patients can stay in the hospital with a probability of 0.7. The total number of days in the Hospital state is expected to be fewer than Sick.
Dead (D): Once patients transition to the Dead state, they remain there indefinitely, so over time, more patients will accumulate in this state, contributing to the overall number of days spent in the Dead state.

The Healthy (H) state should dominate in terms of the total number of days spent.
The Dead (D) state will gradually accumulate as more patients transition into it and stay there.
Sick (S0) and Hospital (SH) will likely show fewer total days compared to the Healthy (H) and Dead (D) states, but they will still play an important role in the patient trajectory.
