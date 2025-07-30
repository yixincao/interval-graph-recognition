package util;

import java.util.Arrays;
import java.security.SecureRandom;
import java.util.stream.IntStream;

/**
 * 
 * @author Yixin Cao and Chenxi Liu (May 2025)
 *
 * Simple utility functions for permutations.
 */
public class Permutation {
    /**
     * To generate a random permutation of [0..n-1] by shuffling the identity permutation.
     * 
     * Starting with the array of the identity permutation, we iterate through the array in a for loop. 
     * For each index i, we swap a[i] with a[r], with a randomly generated index r.
     * 
     * @param n: the number of elements
     * @return a random permutation of [0..n-1]
     */
    public static int[] shuffle(int n) {
        int[] a = IntStream.range(0, n).toArray(); //IntStream.rangeClosed(0, n-1).toArray();
        SecureRandom rand = new SecureRandom();
        for (int i = 0; i < n; i++) {
            int r = rand.nextInt(n);
            int temp = a[r];
            a[r] = a[i];
            a[i] = temp;
        }
        assert (isPermutation(a));
        return a;
    }

    /**
     * Check whether an integer array of length n is a permutation of [0..n-1].
     *
     * We create a boolean array {@code has}, where {@code has[i]} means {@code sigma} contains number i. They are initially all {@code false}.
     * We iterate through {@code sigma} in a for loop.
     * When element i is first met, {@code has[i]} is set to {@code true}.
     * If it is met again, return fasle: {@code sigma} cannot be a permutation.
     *
     * @param sigma: an integer array
     * @return {@code true} if sigma is a permutation of [0..n-1]
     **/
    public static boolean isPermutation(int[] sigma) {
        int n = sigma.length;
        boolean[] has = new boolean[n];
        // Arrays.fill(x, false);
        for (int i : sigma) {
            if (has[i] == true) 
                return false;
            has[i] = true;
        }
         return true;
    }

    /**
     * To reverse the elements within an array.
     * This method is in-place: it modifies the given order.
     * 
     * [1, 0, 3, 4, 2] -> [2, 4, 3, 0, 1]
     * 
     * @param a: the initial permutation
     */
    public static void reverse(int[] a) {
        assert (isPermutation(a));
        int n = a.length;
        for(int i = 0; i < n / 2; i++) {
            int temp = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = temp;
        }
     }

    /**
     * To compute the inverse of a permutation.
     * Note its difference from {@code reverse}.
     *
     * [1, 0, 3, 4, 2] -> [1, 0, 4, 2, 3]
     *
     * inverse[sigma[i]] = sigma[inverse[i]] = i for all i in [0..n-1].
     * 
     * @param sigma: the initial permutation
     * @return the inverse of sigma (\sigma^{-1})
     */
    public static int[] inversePermutation(int[] sigma) {
        assert (isPermutation(sigma));
        int n = sigma.length;
        int[] newlist = new int[n];
        for (int i = 0; i < n; i++) {
            newlist[sigma[i]] = i;
        }
        assert (isPermutation(newlist));
        return newlist;
    }

    /**
     * Random list of M non-negative integers whose sum is N
     * 
     * Given two integers M and N, the task is to create a list of M positive
     * integers whose sum is N.
     * Approach: To get a complete random list of integers, create an array of size
     * M where every element is initialised with 1. Now run a loop from 0 to N - M
     * -1 and increment any randomly chosen element from the array by 1 using the
     * rand() function. This way, sum of the resultant list will be N.
     * 
     * https://www.geeksforgeeks.org/random-list-of-m-non-negative-integers-whose-sum-is-n/
     * 
     * @param m - a list of M numbers
     * @param n - sum of n
     * @return
     */
    public static int[] randomList(int m, int n) {
        // Create an array of size m where
        // every element is initialized to 1
        int arr[] = new int[m];
        Arrays.fill(arr, 1);

        // To make the sum of the final list as n
        for (int i = 0; i < n - m; i++) {

            // Increment any random element
            // from the array by 1
            arr[(int) (Math.random() * m)]++;
        }

        return arr;

    }
}
