public class Main {
    public static void main(String[] args) {
        NeedlemanWunsch n1 = new NeedlemanWunsch("UUAGG", "CGGCC", false);
        n1.alignStrands(); 
    }
}
