package util;

public class Node<T> {
    public T element;
    public Node<T> next;
    public Node<T> previous;

    public Node(T e) {
        element = e;
        next = previous = null;
    }

}
