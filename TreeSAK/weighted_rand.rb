# /usr/bin/env ruby

# all credits to Wolfgang Teuber (https://gitlab.com/knugie)

def weighted_rand(weights = {})
  raise 'Probabilities must sum up to 1' unless weights.values.inject(&:+) == 1.0
  raise 'Probabilities must not be negative' unless weights.values.all? { |p| p >= 0 }
  # Do more sanity checks depending on the amount of trust in the software component using this method,
  # e.g. don't allow duplicates, don't allow non-numeric values, etc.
  
  # Ignore elements with probability 0
  weights = weights.reject { |k, v| v == 0.0 }   # e.g. => {"a"=>0.4, "b"=>0.4, "c"=>0.2}

  # Accumulate probabilities and map them to a value
  u = 0.0
  ranges = weights.map { |v, p| [u += p, v] }   # e.g. => [[0.4, "a"], [0.8, "b"], [1.0, "c"]]

  # Generate a (pseudo-)random floating point number between 0.0(included) and 1.0(excluded)
  u = rand   # e.g. => 0.4651073966724186
  
  # Find the first value that has an accumulated probability greater than the random number u
  ranges.find { |p, v| p > u }.last   # e.g. => "b"
end
