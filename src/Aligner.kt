import java.lang.Math.max

// Class for alignment algorithms
class Aligner (seqj: String,
               seqi: String,
               method: String,
               gap_open: Int,
               gap_extend: Int,
               gap_double: Int,
               max_hits: Int,
               matrix: SubMatrix) {

    // Data structure for results
    data class AlignResults (val seq1: String, val seq2: String,
                             var start1: Int, var start2: Int,
                             val end1: Int, val end2: Int, val score: Int,
                             var n_gaps1: Int, var n_gaps2: Int, var n_mismatches: Int) {
        override fun toString(): String {
            return seq1 + "\n" + seq2
        }
    }

    // Early initialised variables
    val METHODS = arrayOf("global", "local", "glocal")
    val GAP_CHAR = '-'
    val NONE = 0
    val LEFT = 1
    val UP = 2
    val DIAG = 3

    // Late initialised variables
    val max_i: Int
    val max_j: Int
    val seq_i: String
    val seq_j: String
    val flip: Boolean

    // Results
    val f_score: Array<IntArray>
    val pointer: Array<IntArray>
    val results: List<AlignResults>


    private fun assignScores(meth: String, gapExt: Int, gapOpen: Int, gapDoub: Int, matr: SubMatrix): Unit {
        // Initial matrix
        val i_score = Array(max_i + 1) { IntArray(max_j + 1) { Int.MIN_VALUE } }
        val j_score = Array(max_i + 1) { IntArray(max_j + 1) { Int.MIN_VALUE } }
        when (meth) {
            "global" -> {
                for (i in 1 until max_j + 1) {pointer[0][i] = LEFT}
                for (i in 1 until max_i + 1) {pointer[i][0] = UP}
                for (i in 1 until max_j + 1) {f_score[0][i] = gapOpen + gapExt * (i - 1)}
                for (i in 1 until max_i + 1) {f_score[i][0] = gapOpen + gapExt * (i - 1)}
            } "glocal" -> {
                for (i in 1 until max_j + 1) {pointer[0][i] = LEFT}
                for (i in 1 until max_j + 1) {f_score[0][i] = gapOpen + gapExt * (i - 1)}
            }
        }

        // Fill matrix
        for (i in 1 until max_i + 1) {
            val charI: Char = seq_i[i - 1]
            for (j in 1 until max_j + 1) {
                val charJ = seq_j[j - 1]

                // Assign I, J, F matrix at index (i, j)
                i_score[i][j] = listOf(
                        f_score[i][j - 1] + gapOpen,
                        if (i_score[i][j - 1] == Int.MIN_VALUE) Int.MIN_VALUE else i_score[i][j - 1] + gapExt,
                        if (j_score[i][j - 1] == Int.MIN_VALUE) Int.MIN_VALUE else j_score[i][j - 1] + gapDoub
                ).max() ?: 0
                j_score[i][j] = listOf(
                        f_score[i - 1][j] + gapOpen,
                        if (i_score[i - 1][j] == Int.MIN_VALUE) Int.MIN_VALUE else i_score[i - 1][j] + gapExt,
                        if (j_score[i - 1][j] == Int.MIN_VALUE) Int.MIN_VALUE else j_score[i - 1][j] + gapDoub
                ).max() ?: 0
                val diag_score = f_score[i - 1][j - 1] + matr.getScore(charJ, charI)
                val left_score = i_score[i][j]
                val up_score = j_score[i][j]
                val max_score = listOf(diag_score, up_score, left_score).max() ?: 0
                f_score[i][j] = if (meth == "local") max(0, max_score) else max_score

                // Assign pointer matrix according to alignment method
                when (meth) {
                    in listOf("local", "glocal") -> {
                        when (max_score) {
                            up_score -> pointer[i][j] = UP
                            left_score -> pointer[i][j] = LEFT
                            diag_score -> pointer[i][j] = DIAG
                        }
                    } "global" -> {
                        when (max_score) {
                            up_score -> pointer[i][j] = UP
                            left_score -> pointer[i][j] = LEFT
                            else -> pointer[i][j] = DIAG
                        }
                    }
                }
            }
        }
    }


    private fun getPairs(limitVal: Int, meth: String) : List<Pair<Int, Int>>{
        val outList = mutableListOf<Pair<Int, Int>>()

        // Find indexes
        when (meth) {
            "global" -> outList.add(Pair(max_i, max_j))
            "local" -> {
                val maxVals: List<Int> = f_score.map{ it -> it.asList().max() ?: 0}
                for (i in 1 until max_i + 1) {
                    for (j in 1 until max_j) { if (maxVals.max() == f_score[i][j]) outList.add(Pair(i, j)) }
                }
            }
            "glocal" -> {
                val endVals: List<Int> = f_score.map{ it -> it.last() }
                for (i in 1 until max_i + 1) {
                    if (endVals.max() == f_score[i].last()) outList.add(Pair(i, max_j))
                }
            }
        }

        // Return operation
        return outList.shuffled().take(limitVal)
    }


    private fun getResults(pairs: List<Pair<Int, Int>>, meth: String) : List<AlignResults> {
        val resultArray = mutableListOf<AlignResults>()
        for (pair in pairs) {
            // Get indexes
            var i = pair.first
            var j = pair.second

            // Set constant variables
            val score = f_score[i][j]
            val align_j = mutableListOf<Char>()
            val align_i = mutableListOf<Char>()
            val end_i: Int
            val end_j: Int
            when (meth) {
                "global" -> {
                    end_i = max_i
                    end_j = max_j
                } else -> {
                    end_i = i
                    end_j = j
                }
            }

            // Set mutable variables
            var p = pointer[i][j]
            var n_gaps_i = 0
            var n_gaps_j = 0
            var n_mmatch = 0

            // Path traversal algorithm with the pointer matrix
            while (p != NONE) {
                when (p) {
                    DIAG -> {
                        i -= 1
                        j -= 1
                        val ichar = seq_i[i]
                        val jchar = seq_j[j]
                        if (ichar != jchar) { n_mmatch += 1 }
                        align_j.add(jchar)
                        align_i.add(ichar)
                    } LEFT -> {
                        j -= 1
                        align_j.add(seq_j[j])
                        if ( align_i.isEmpty() || align_i.last() != GAP_CHAR ) { n_gaps_i += 1 }
                        align_i.add(GAP_CHAR)
                    } UP -> {
                        i -= 1
                        align_i.add(seq_i[i])
                        if ( align_j.isEmpty() || align_j.last() != GAP_CHAR ) { n_gaps_j += 1 }
                        align_j.add(GAP_CHAR)
                    } else -> Exception("Invalid pointer")
                }
                p = pointer[i][j]
            }

            // Append result to results list
            val str_align_i = align_i.reversed().joinToString(separator = "") { it.toString() }
            val str_align_j = align_j.reversed().joinToString(separator = "") { it.toString() }
            val thisResult: AlignResults
            if (flip) {
                thisResult = AlignResults(str_align_i, str_align_j, i, j, end_i, end_j,
                        n_gaps_i, n_gaps_j, n_mmatch, score)
            } else {
                thisResult = AlignResults(str_align_j, str_align_i, j, i, end_j, end_i,
                        n_gaps_j, n_gaps_i, n_mmatch, score)
            }
            resultArray.add(thisResult)
        }
        return resultArray
    }


    init {
        // 0. Assertions
        assert(max_hits > 0)
        assert(method in METHODS)

        // 1. Initialise variables
        if (seqj.length > seqi.length){
            flip = true
            seq_i = seqj
            seq_j = seqi
            max_i = seqj.length
            max_j = seqi.length
        } else {
            flip = false
            seq_i = seqi
            seq_j = seqj
            max_i = seqi.length
            max_j = seqj.length
        }
        f_score = Array(max_i + 1) { IntArray(max_j + 1) {0} }
        pointer = Array(max_i + 1) { IntArray(max_j + 1) {NONE} }

        // 2. Fill f-score matrix and pointer matrix
        assignScores(method, gap_extend, gap_open, gap_double, matrix)

        // 3. Find traceback coordinates
        val ij_pairs = getPairs(max_hits, method)

        // 4. Find paths across f-score matrix
        results = getResults(ij_pairs, method)
    }
}