import unittest
import eval_methods
from eval_methods import Seq_has


class Test_is_spurious_by_max_repeat_len_and_min_distance_seq(unittest.TestCase):
    def test_is_spurious_by_max_repeat_len_and_min_distance_seq_has_no_neighbour(self):

        def f(s, max_repeat_len=3, min_distance=6):
            return eval_methods.is_spurious_by_max_repeat_len_and_min_distance(s, max_repeat_len, min_distance,
                                                                               Seq_has.NO_NEIGHBOUR)
        self.assertTrue(f("(CAG)_2 NNNNNN(CAG)_2 "))
        self.assertFalse(f("(CAG)_2 NNNNNN(CAG)_4 "))
        self.assertFalse(f("(CAG)_2 NNNNN(CAG)_2 "))

    def test_is_spurious_by_max_repeat_len_and_min_distance_seq_has_left_neighbour(self):

        def f(s, max_repeat_len=3, min_distance=6):
            return eval_methods.is_spurious_by_max_repeat_len_and_min_distance(s, max_repeat_len, min_distance,
                                                                               Seq_has.LEFT_NEIGHBOUR)
        self.assertTrue(f("NNNNNN(CAG)_2 N"))
        self.assertFalse(f("NNNNN(CAG)_2 "))
        self.assertFalse(f("NNNNNN(CAG)_4 NNNNNN"))

    def test_is_spurious_by_max_repeat_len_and_min_distance_seq_has_right_neighbour(self):

        def f(s, max_repeat_len=3, min_distance=6):
            return eval_methods.is_spurious_by_max_repeat_len_and_min_distance(s, max_repeat_len, min_distance,
                                                                               Seq_has.RIGHT_NEIGHBOUR)
        self.assertTrue(f("N(CAG)_2 NNNNNN"))
        self.assertFalse(f("(CAG)_2 NNNNN"))
        self.assertFalse(f("NNNNNN(CAG)_4 NNNNNN"))

    def test_is_spurious_by_max_repeat_len_and_min_distance_seq_has_both_neighbours(self):

        def f(s, max_repeat_len=3, min_distance=6):
            return eval_methods.is_spurious_by_max_repeat_len_and_min_distance(s, max_repeat_len, min_distance,
                                                                               Seq_has.BOTH_NEIGHBOURS)
        self.assertTrue(f("NNNNNN(CAG)_2 NNNNNN"))
        self.assertFalse(f("NNNNNN(CAG)_2 NNNNN"))
        self.assertFalse(f("NNNNN(CAG)_2 NNNNNN"))
        self.assertFalse(f("NNNNNN(CAG)_4 NNNNNN"))


class Test_seq_conforms_with_category(unittest.TestCase):

    categories = [["(CAG)", "TAG", "(CAG)"],
                  ["(CAG)", "CAA", "(CAG)"],
                  ["(CAG)", "CCG", "(CAG)"],
                  ["(CAG)"]]
    max_repeat_len = 3
    min_distance = 6

    def test_seq_conforms_with_category(self):

        def f(x):
            col = {"no_flanks": x}
            return list(eval_methods.seq_conforms_with_category(col, Test_seq_conforms_with_category.categories,
                                                                lambda seq, neighbour:
                                                                eval_methods.is_spurious_by_max_repeat_len_and_min_distance(
                                                                    seq,
                                                                    Test_seq_conforms_with_category.max_repeat_len,
                                                                    Test_seq_conforms_with_category.min_distance,
                                                                    neighbour)))

        self.assertListEqual(f("(CAG)_2 TAG(CAG)_12 "), [1, 0, 0, 1, True])
        self.assertListEqual(f("(CAG)_2 TAG(CAG)_12 CAA(CAG)_2 "), [1, 1, 0, 1, True])
        self.assertListEqual(f("(AAA)_2 NNNNNN(AAA)_2 "), [0, 0, 0, 0, True])
        self.assertListEqual(f("(AAA)_4 NNNNNN(AAA)_2 "), [0, 0, 0, 0, False])
        self.assertListEqual(f("(CAG)_2 TAG(CAG)_12 CAA(CAG)_2 NNNNN(TTT)_3 "), [1, 1, 0, 1, False])
        self.assertListEqual(f("(CAG)_2 TAG(CAG)_12 CAA(CAG)_2 NNNNNN(TTT)_3 "), [1, 1, 0, 1, True])
        self.assertListEqual(f("(TTT)_3 NNNNN(CAG)_2 TAG(CAG)_12 CAA(CAG)_2 "), [1, 1, 0, 1, False])
        self.assertListEqual(f("(TTT)_3 NNNNNN(CAG)_2 TAG(CAG)_12 CAA(CAG)_2 "), [1, 1, 0, 1, True])
        self.assertListEqual(f("(CAG)_2 TAG(CAG)_12 NNNNNN(AAA)_3 NNNNNN(CAG)_2 CAA(CAG)_2 NNNNNN(TTT)_3 "), [1, 1, 0, 1, True])
        self.assertListEqual(f("(CAG)_2 TAG(CAG)_12 NNNNNN(AAA)_3 NNNNN(CAG)_2 CAA(CAG)_2 NNNNNN(TTT)_3 "), [1, 1, 0, 1, False])
        self.assertListEqual(f("(CAG)_2 TAG(CAG)_12 NNNNN(AAA)_3 NNNNNN(CAG)_2 CAA(CAG)_2 NNNNNN(TTT)_3 "), [1, 1, 0, 1, False])
        self.assertListEqual(f("(CAG)_2 TAG(CAG)_12 NNNNNN(AAA)_4 NNNNNN(CAG)_2 CAA(CAG)_2 NNNNNN(TTT)_3 "), [1, 1, 0, 1, False])
        self.assertListEqual(f("(CCA)_2 NNNNNN(CAC)_2 NNNNNN(CAG)_15 NNNNNN(CCT)_3 "), [0, 0, 0, 1, True])
        self.assertListEqual(f("(CCA)_2 NNNNNN(CAC)_2 NNNNN(CAG)_15 NNNNNN(CCT)_3 "), [0, 0, 0, 1, False])
        self.assertListEqual(f("(CCA)_2 NNNNNN(CAG)_3 TAG(CAG)_4 NNNNN(CAG)_15 NNNNNN(CCT)_3 "), [1, 0, 0, 1, False])
        self.assertListEqual(f("(CCA)_2 NNNNNN(CAG)_3 TAG(CAG)_4 NNNNNN(CAG)_15 NNNNNN(CCT)_3 "), [1, 0, 0, 1, True])
        self.assertListEqual(f("(CCA)_2 NNNNNN(CAG)_3 TAG(CAG)_4 NNNNNN(CAG)_15 NNNNNN(CAG)_6 NNNN"), [1, 0, 0, 1, True])


if __name__ == '__main__':
    unittest.main()
