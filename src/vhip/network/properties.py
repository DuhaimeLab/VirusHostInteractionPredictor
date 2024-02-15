"""Adjacency matrix class that takes for input a numpy 2d array."""

from typing import List

import numpy as np
from numpy.typing import NDArray


class AdjacencyMatrix:
    """Define class for virus-host interaction networks (VHINs).

    Args:
        matrix (np.array): Array representing interactions.
    """

    def __init__(self, matrix: NDArray[np.int_]) -> None:
        """Initializes class variables."""
        self.adj = matrix

        self.sorted = False
        self.shape = matrix.shape

    def nestedness(self) -> None:
        """Calculate nestedness using the NODF algorithm for the array."""
        # sort matrix if it has not been done
        if self.sorted is False:
            self.sort()

        N_row = 0
        N_col = 0

        # generate list of rows to compare
        rows_to_compare = self.pairs(axis=0)
        for x, y in rows_to_compare:
            pair1 = self.adj[x,]
            pair2 = self.adj[y,]

            N_row += self.compare(pair1, pair2)

        # generate list of columns to compare
        cols_to_compare = self.pairs(axis=1)
        for x, y in cols_to_compare:
            pair1 = self.adj[:, x]
            pair2 = self.adj[:, y]
            print(pair1, pair2)

            N_col += self.compare(pair1, pair2)  # pyright: ignore

        print(N_row, N_col)
        self.nodf = (N_row + N_col) / 2

    def sort(self) -> None:
        """Sort matrix."""
        # calculate sums of 1s for rows and columns, and stores as a list of tuples
        # first element of tuple: sum | second element: index
        self.sum_rows: List[tuple[int, int]] = []
        self.sum_cols: List[tuple[int, int]] = []

        for i, row in enumerate(self.adj):
            self.sum_rows.append((sum(row), i))
        for i, col in enumerate(self.adj.T):
            self.sum_cols.append((sum(col), i))

        # descending order for sum_rows and sum_cols
        self.sum_rows = sorted(self.sum_rows, reverse=True)
        self.sum_cols = sorted(self.sum_cols, reverse=True)

        # unpack second element of list of tuples to get desired order of indices
        new_row_order = [x[1] for x in self.sum_rows]
        new_col_order = [x[1] for x in self.sum_cols]

        # rearrange matrix
        self.adj = self.adj[new_row_order,]
        self.adj = self.adj[:, new_col_order]
        self.sorted = True

    def pairs(self, axis: int = 0) -> list[tuple[int, int]]:
        """Determine all possible i-j pairs.

        Args:
            axis (int): Axis to be used when determining all pairs.
        """
        lst: List[tuple[int, int]] = []
        for i in range(0, self.shape[axis]):
            for j in range(i + 1, self.shape[axis]):
                lst.append((i, j))
        return lst

    def compare(self, x: list[int], y: list[int]) -> float:
        """Compare two lists containing 0 and 1.

        Args:
            x (list[int]): first list
            y (list[int]): second list
        """
        if sum(x) <= sum(y):
            val = 0
        else:
            counter = 0
            total = 0
            for i, j in zip(x, y):
                if i == 1 and j == 1:
                    counter += 1
                    total += 1
                elif i == 0 and j == 1:
                    total += 1
            val = counter / total

        return val * 10


# TODO: Export code below to be its own test

unsorted = np.array(
    [
        [0, 0, 0, 1, 1],
        [0, 1, 1, 1, 0],
        [0, 0, 0, 1, 1],
        [0, 0, 1, 1, 1],
        [1, 1, 1, 0, 1],
    ]
)

test = AdjacencyMatrix(unsorted)
test.nestedness()
print(test.nodf)
