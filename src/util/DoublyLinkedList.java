package util;

/**
 * A simple implementation of doubly linked list data structure.
 */

public class DoublyLinkedList<T> {
    public Node<T> head;
    public Node<T> tail;
    // when was the linked list created
    public int time;
    private int size;

    public DoublyLinkedList() {
        this(-1);
    }

    public DoublyLinkedList(int time) {
        this.head = null;
        this.tail = null;
        size = 0;
        this.time = time;
    }

    public boolean isEmpty() {
        return size == 0;
    }

    public void err() {
        System.out.println("Oops...");
    }

    public void insertAtBeginning(Node<T> temp) {
        if (head == null) {
            head = temp;
            tail = temp;
        } else {
            temp.next = head;
            head.previous = temp;
            head = temp;
        }
        size++;
    }

    public void insertBefore(Node<T> old, Node<T> e) {
        if (this.head == old) {
            insertAtBeginning(e);
        } else {

            e.next = old;
            e.previous = old.previous;
            old.previous.next = e;
            old.previous = e;

            size++;
        }
    }

    public void insertAtEnd(Node<T> e) {
        if (tail == null) {
            head = e;
            tail = e;
        } else {
            tail.next = e;
            e.previous = tail;
            tail = e;
        }
        size++;
    }

    public void deleteAtBeginning() {
        if (head == null) {
            return;
        }

        if (head == tail) {
            head = null;
            tail = null;
            size--;
            return;
        }

        Node<T> temp = head;
        head = head.next;
        head.previous = null;
        temp.next = null;
        size--;
    }

    public void delete(Node<T> node) {
        if (head == null) {
            return;
        }

        if (this.head == node) {
            deleteAtBeginning();
            return;
        }

        if (this.tail == node) {
            deleteAtEnd();
            return;
        }

        node.previous.next = node.next;
        node.next.previous = node.previous;
        node.previous = null;
        node.next = null;
        size--;
    }

    public void deleteAtEnd() {
        if (tail == null) {
            return;
        }

        if (head == tail) {
            head = null;
            tail = null;
            size--;
            return;
        }

        Node<T> temp = tail;
        tail = tail.previous;
        tail.next = null;
        temp.previous = null;
        size--;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(("["));
        Node<T> temp = this.head;
        while (temp != null) {
            sb.append(temp.element.toString());
            temp = temp.next;
            if (temp != null)
                sb.append(" -> ");
        }
        return sb.append("]").toString();
    }
    
    public void display() {
        Node<T> temp = this.head;
        while (temp != null) {
            System.out.print(temp.element.toString() + " --> ");
            temp = temp.next;
        }
        System.out.println("NULL");
    }

}
