fun main(args : Array<String>) {
    // Arguments
    val seqj = "ACTG"
    val seqi = "CTTTTCGCGG"
    val method = "glocal"
    val gap_open = -5
    val gap_extend = -5
    val gap_double = -5
    val max_hits = 3

    // Substitution matrix
    var matrix = SubMatrix()

    // Alignment
    val aligned = Aligner(seqj, seqi, method, gap_open, gap_extend, gap_double, max_hits, matrix)
    for (l in aligned.results){ println(l) }
}