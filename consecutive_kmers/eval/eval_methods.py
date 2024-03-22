import re
from typing import List, Callable
from enum import Enum
import pandas as pd
import itertools


class Seq_has(Enum):
    NO_NEIGHBOUR = 0
    LEFT_NEIGHBOUR = 1
    RIGHT_NEIGHBOUR = 2
    BOTH_NEIGHBOURS = 3


def seq_conforms_with_category(row,
                               categories: List[List[str]],
                               is_spurious: Callable[[str, Seq_has], bool] = lambda x: False
                               ) -> pd.Series:
    """
    Parameters:
    row: a row from the dataframe
    categories: a list of lists of strings, each list of strings represents a category of tandem repeats
    is_spurious: a function that takes a string and returns a boolean, if repeats in a str are spurious

    Returns: a series of booleans, each representing if the row's repeat_representation conforms with a category. The
    last element is a boolean whether the seq is completely categorized or not
    """
    seq = row["no_flanks"]  # if this is changed then i must adapt middles_are_spurious
    found_categories = [0] * len(categories)
    positions_with_found_categories = [0] * len(seq)
    for i, category in enumerate(categories):
        # find all start and stops, see if all tandem repeats in the string are covered
        pattern = (p if p[0] != "(" else p.replace("(", "\(").replace(")", "\)") + "\_\d+\s" for p in category)
        pattern = "(?=(" + "".join(pattern) + "))"  # to find overlapping matches
        matches = re.finditer(pattern, seq)
        for hit in matches:

            start = hit.start()
            end = hit.end() + len(hit.group(1))

            found_categories[i] = 1
            positions_with_found_categories[start:end] = [1] * (end - start)

    mapped_string = "".join(str(x) for x in positions_with_found_categories)
    # print("mapped string ", mapped_string)
    matches = [match for match in re.finditer("0+", mapped_string)]
    if len(matches) == 0:
        return pd.Series(found_categories + [True])  # seq is complete categorized
    found_not_spurious = False
    for match in matches:  # check all un-categorized parts
        if match.start() == 0 and match.end() == len(seq):
            found_not_spurious = found_not_spurious or not is_spurious(seq, Seq_has.NO_NEIGHBOUR)
        elif match.start() == 0:
            found_not_spurious = found_not_spurious or not is_spurious(seq[0:match.end()], Seq_has.RIGHT_NEIGHBOUR)
        elif match.end() == len(seq):
            found_not_spurious = found_not_spurious or not is_spurious(seq[match.start():], Seq_has.LEFT_NEIGHBOUR)
        else:
            found_not_spurious = found_not_spurious or not is_spurious(
                seq[match.start():match.end()], Seq_has.BOTH_NEIGHBOURS)

    return pd.Series(found_categories + [not found_not_spurious])


def is_spurious_by_max_repeat_len_and_min_distance(seq,
                                                   max_repeat_len: int,
                                                   min_distance: int,
                                                   relation_to_neighbour: Seq_has,
                                                   ) -> bool:
    """
    Parameters:
    str: the string to check
    max_repeat_len: if repeat is shorter than this, it is considered spurious
    min_distance: if the distance between two repeats is greater than this, it is considered spurious
    relation_to_neighbour: if Seq_has.NO_NEIGHBOUR: then the seq is considered to have no neighbouring sequences
                           if Seq_has.LEFT_NEIGHBOUR: then the seq is considered to have a neighbour to the left that ends in a repeat
                           if Seq_has.RIGHT_NEIGHBOUR: then the seq is considered to have a neighbour to the right that starts with a repeat
                           if Seq_has.BOTH_NEIGHBOURS: then the seq is considered to have both neighbours
    """
    # print(seq)
    for x in re.finditer("\(\w+\)_(\d+)\s", seq):
        if int(x.group(1)) > max_repeat_len:
            # print("repeat is too long")
            return False

    if (x := re.search("\(\w+\)_\d+\s(\w*)\(", seq)):  # finding to repeats which are close to each other in the not yet categorized string
        if len(x.group(1)) < min_distance:
            # print("repeats are to close")
            return False

    if relation_to_neighbour == Seq_has.LEFT_NEIGHBOUR or relation_to_neighbour == Seq_has.BOTH_NEIGHBOURS:
        if (x := re.match("(\w*)\(", seq)):
            if len(x.group(1)) < min_distance:
                # print("Left neighbour is too close")
                return False
    if relation_to_neighbour == Seq_has.RIGHT_NEIGHBOUR or relation_to_neighbour == Seq_has.BOTH_NEIGHBOURS:
        if (x := re.search("\s(\w*)$", seq)):
            if len(x.group(1)) < min_distance:
                # print("Right neighbour is too close")
                return False

    if relation_to_neighbour == Seq_has.BOTH_NEIGHBOURS:
        if len(seq) < min_distance:
            # print("seq is too short")
            return False
    return True


class Color_print_triplets:

    background_colors = [40, 41, 42, 43, 44, 45, 46, 47]
    font_colors = [30, 31, 32, 33, 34, 35, 36, 37]

    def __init__(self):
        color_combinations = set(itertools.product(self.background_colors, self.font_colors))
        self.colors_for_triplets = {"GAC": (41, 30),
                                    "TAG": (42, 30)}
        for triplet, color_combination in self.colors_for_triplets.items():
            color_combinations.remove(color_combination)
        for triplet in itertools.product("ACGT", repeat=3):
            triplet = "".join(triplet)
            if triplet not in self.colors_for_triplets.keys():
                self.colors_for_triplets[triplet] = color_combinations.pop()

    def color_triplets(self, seq, expand=False):
        pattern = r'\(([A-Z]+)\)_(\d+)\s'

        def repl(m):
            triplet = m.group(1)
            new_repeat = m.group(1) if not expand else triplet * int(m.group(2))
            colorized_triplet = f"\033[1;{self.colors_for_triplets[triplet][1]};{self.colors_for_triplets[triplet][0]}m{new_repeat}\033[0m"
            return colorized_triplet
        return re.sub(pattern, repl, seq)


class Color_pd:
    # this is currently not used
    def color_negative_red(val):
        if val < 1:
            font_color = 'white'
            background_color = 'red'
        else:
            font_color = 'black'
            background_color = 'white'
        return 'color: %s; background-color: %s' % (font_color, background_color)

    def highlight_chars_with_A(text):
        styled_text = ''
        for char in text:
            if char == 'A':
                styled_text += '<span style="color: white; background-color: red;">%s</span>' % char
            else:
                styled_text += char
        return styled_text

    # df.style.applymap(color_negative_red, subset=["frame"])
    # df['repeat_representation'] = df['repeat_representation'].apply(highlight_chars_with_A)
