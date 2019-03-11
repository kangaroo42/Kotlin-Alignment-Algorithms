import java.lang.Math.max
import java.lang.Math.min

// Class for scoring matrix

class SubMatrix {

    val residues: String = "ACGT"

    private val scores: Array<IntArray> = arrayOf(
            /* A */ intArrayOf(4),
            /* C */ intArrayOf(0, 9),
            /* G */ intArrayOf(0, -3, 6),
            /* T */ intArrayOf(0, -1, -2, 5)
    )

    fun getScore(charA: Char, charB: Char): Int {
        val indexA = residues.indexOf(charA)
        val indexB = residues.indexOf(charB)
        val minIndex = min(indexA, indexB)
        val maxIndex = max(indexA, indexB)
        return scores[maxIndex][minIndex]
    }

    // Check matrix validity
    init {
        assert(residues.length == scores.size)
        for (i in 0 until residues.length) {
            assert(scores[i].size == i + 1)
        }
    }
}