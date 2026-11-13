import pandas as pd
from collections import Counter



def BWT_transformer(input_string:str) -> pd.DataFrame:
    current_string = input_string + '$'
    current_string_len = len(current_string)
    rotations = list()
    rotations.append(current_string)
    for i in range(current_string_len-1):
        current_string = current_string[-1]+current_string[:-1]
        rotations.append(current_string)

    rotations_dataframe = pd.DataFrame({'rotations':rotations,'sorted_rotations':sorted(rotations)})

    F_BWT = [string[0] for string in sorted(rotations)]
    L_BWT = [string[-1] for string in sorted(rotations)]

    count_of_F_BWT = Counter(F_BWT)
    print(count_of_F_BWT)
    F_BWT_counted = []

    for i in sorted(F_BWT):
        F_BWT_counted.append(i + str(count_of_F_BWT[i]))
        count_of_F_BWT[i] -= 1
    F_BWT_counted = sorted(F_BWT_counted)
    print(F_BWT_counted)


    count_of_L_BWT = Counter(L_BWT)
    L_BWT_counted = []

    for i in L_BWT[::-1]:
        L_BWT_counted.append(i+str(count_of_L_BWT[i]))
        count_of_L_BWT[i] -=1
    L_BWT_counted = L_BWT_counted[::-1]
    print(L_BWT_counted)

    BWT_dataframe = pd.DataFrame({'first' : F_BWT_counted,'last' : L_BWT_counted})
    print(BWT_dataframe)
    return(BWT_dataframe)


def BWT_inverting(BWT_dataframe:pd.DataFrame) -> str:
    L2F_mapping = {k:v for k, v in zip(BWT_dataframe['first'].values,range(len(BWT_dataframe)))}
    # print(L2F_mapping)
    BWT_dataframe['L2F'] = BWT_dataframe.apply(lambda row : L2F_mapping[row['last']],axis=1)
    print(BWT_dataframe)

def BWT_LF_mapping (BWT_dataframe : pd.DataFrame, pattern : str) -> list:

    F_L_dict = dict(zip(BWT_dataframe['first'], BWT_dataframe['last']))
    F_dict = {label: idx for idx, label in enumerate(BWT_dataframe['first'])}
    L_dict = {label: idx for idx, label in enumerate(BWT_dataframe['last'])}

    current_labels = set(BWT_dataframe['last'])

    for letter in reversed(pattern):
        matched = {label for label in current_labels if label[0] == letter}
        if not matched:
            return []
        current_labels = {F_L_dict[label] for label in matched}
        print(current_labels)

    position = []

    for label in current_labels:
        row = L_dict[label]
        step = 0
        while BWT_dataframe.iloc[row]['first'][0] != '$':
            row = F_dict[BWT_dataframe.iloc[row]['last']]
            step += 1
        position.append(step-1)

    print(position)
    return sorted(position)



BWT_test_dataframe = BWT_transformer('BANANA')
BWT_LF_mapping(BWT_test_dataframe,'A')

#BWT_inverting(BWT_test_dataframe)

