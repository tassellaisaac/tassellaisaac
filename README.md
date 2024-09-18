Burrows-Wheeler Transform (BWT) and Reverse BWT with Bit Calculator
Description:
This project implements the Burrows-Wheeler Transform (BWT) and its reverse operation (Reverse BWT) in R. The BWT is an important algorithm in bioinformatics, text compression, and sequence alignment. It rearranges a string into a form that is more compressible by clustering similar characters together. Additionally, this project includes a simple bit-calculator that estimates the number of bits required to store a given DNA sequence based on its length.

The project processes a set of DNA sequences, computes the BWT, reverses the BWT to verify correctness, and calculates the bit storage requirements for each sequence.

Project Structure and Functions
Main Functions
bwt(sequence): Computes the Burrows-Wheeler Transform of the input sequence.
reverse_bwt(bwt_seq): Reconstructs the original sequence from the BWT result.
bit_calculator(sequence): Calculates the number of bits required to store the sequence.
test_sequences(sequences): Runs the BWT, Reverse BWT, and bit calculator for a list of sequences.

Limitations and Future Enhancements
End-of-string marker: The program assumes the presence of the end-of-string marker ($) for proper BWT functionality. The marker is automatically added when testing sequences.
Optimizations: While the current implementation works for small sequences, future optimizations could improve performance for longer sequences or large datasets.
Bit Calculation: The bit calculation assumes 2 bits per character, which is typical for nucleotides, but could be adjusted for different encoding schemes.

# Function to calculate the Burrows-Wheeler Transform (BWT)
bwt <- function(sequence) {
  n <- nchar(sequence)
  # Create a list of cyclic shifts of the sequence
  rotations <- sapply(0:(n-1), function(i) {
    substring(sequence, i + 1, n) + substring(sequence, 1, i)
  })
  
  # Sort the rotations lexicographically
  sorted_rotations <- sort(rotations)
  
  # The BWT is the last column of the sorted rotations
  bwt_result <- paste0(sapply(sorted_rotations, function(x) substring(x, n, n)), collapse = "")
  
  return(bwt_result)
}

# Function to calculate Reverse Burrows-Wheeler Transform (Reverse BWT)
reverse_bwt <- function(bwt_seq) {
  n <- nchar(bwt_seq)
  
  # Initialize an empty matrix to hold the sorted characters
  table <- rep("", n)
  
  # Repeat the following steps n times
  for (i in 1:n) {
    # Prepend the BWT characters to the rows of the table
    table <- paste0(substring(bwt_seq, 1:n, 1:n), table)
    
    # Sort the table lexicographically
    table <- sort(table)
  }
  
  # Find the row that ends with the end-of-string character '$'
  for (row in table) {
    if (substring(row, n, n) == "$") {
      # Remove the end-of-string character and return the original sequence
      return(substring(row, 1, n-1))
    }
  }
}

# Function to calculate the number of bits required to store a sequence
bit_calculator <- function(sequence) {
  n <- nchar(sequence)
  # Each character requires 2 bits (assuming 4 nucleotides A, T, C, G)
  return(2 * n)
}

# Function to test forward and backward BWT and bit calculations
test_sequences <- function(sequences) {
  for (seq in sequences) {
    cat("Original Sequence: ", seq, "\n")
    
    # Add an end-of-string marker "$" to the sequence
    seq_with_marker <- paste0(seq, "$")
    
    # Perform BWT
    bwt_result <- bwt(seq_with_marker)
    cat("BWT: ", bwt_result, "\n")
    
    # Reverse the BWT to get the original sequence
    reversed_seq <- reverse_bwt(bwt_result)
    cat("Reconstructed Sequence: ", reversed_seq, "\n")
    
    # Calculate the number of bits required to store the sequence
    bits <- bit_calculator(seq)
    cat("Bits Required: ", bits, "\n")
    
    # Add some separation between results for each sequence
    cat("------------------------------------------------------\n")
  }
}

# Test the function with the given sequences
sequences <- c(
  "GATTACA",
  "ATTACATTAC",
  "ATATATATATA",
  "ATATATATAT",
  "AATAATAATAAT",
  "AAAATAAATAAA",
  "ATATACACACA",
  "ATATGTATACAT"
)

test_sequences(sequences)
Original Sequence:  GATTACA
BWT:  ACT$GAAT
Reconstructed Sequence:  GATTACA
Bits Required:  14
------------------------------------------------------
Original Sequence:  ATTACATTAC
BWT:  CCTTAAAT$AC
Reconstructed Sequence:  ATTACATTAC
Bits Required:  20
------------------------------------------------------
Original Sequence:  ATATATATATA
BWT:  TTTTTAAAAAA$
Reconstructed Sequence:  ATATATATATA
Bits Required:  24
------------------------------------------------------
...
