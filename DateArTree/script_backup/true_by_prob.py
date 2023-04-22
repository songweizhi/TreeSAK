import random


def true_by_prob(true_prob):

    true_num = round(1000*true_prob)
    false_num = round(1000*(1 - true_prob))
    option_list = true_num*[True] + false_num*[False]
    decision_to_return = random.choice(option_list)

    return decision_to_return


true_ratio = 0.097
decision_made = true_by_prob(true_ratio)
print(decision_made)

