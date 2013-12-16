#include "HuffWT.h"
#include <queue>
#include <vector>

HuffWT::HuffWT(uchar *s, ulong n, TCodeEntry *codetable, unsigned level) 
    :bitrank(0), left(0), right(0), codetable(0), ch(0), leaf(0)
{
    ch = s[0];
    leaf = false;
    this->codetable = codetable;
    
    bool *B = new bool[n];
    ulong sum=0,i;
    /*
    for (i=0; i< n; i++) {
       printf("%c:", (char)((int)s[i]-128));
       for (r=0;r<codetable[(int)s[i]].bits;r++)
          if (codetable[(int)s[i]].code & (1u <<r))
	     printf("1");
	  else printf("0");
       printf("\n");	  
    }
    printf("\n");
    if (level > 100) return;*/
    for (i=0;i<n;i++) 
       if (codetable[(int)s[i]].code & (1u << level)) {
          B[i] = true;
          sum++;
       }
       else B[i] = false;
    if (sum==0 || sum==n) {
        delete [] B;
	leaf = true;
        return;
    } 
    uchar *sfirst, *ssecond;
    //if (n-sum > 0)  
        sfirst = new uchar[n-sum];
    //if (sum > 0) 
        ssecond = new uchar[sum];
    unsigned j=0,k=0;
   for (i=0;i<n;i++)
        if (B[i]) ssecond[k++] = s[i];
        else sfirst[j++] = s[i];
    ulong *Binbits = new ulong[n/W+1];
    for (i = 0; i < n/W+1; ++i)
        Binbits[i] = 0;
    for (i=0;i<n;i++)
        Tools::SetField(Binbits,1,i,B[i]); 
    delete [] B;
    bitrank = new BitRank(Binbits,n,true);
    //if (j > 0) { 
        left = new HuffWT(sfirst,j,codetable,level+1); 
        delete [] sfirst;
    //}
    //if (k>0) {
        right = new HuffWT(ssecond,k,codetable,level+1); 
        delete [] ssecond;
    //}
}

HuffWT::HuffWT(std::FILE *file, TCodeEntry *ct)
    :bitrank(0), left(0), right(0), codetable(ct), ch(0), leaf(0)
{
    if (std::fread(&leaf, sizeof(bool), 1, file) != 1)
        throw std::runtime_error("HuffWT: file read error (Rs).");
    if (std::fread(&ch, sizeof(uchar), 1, file) != 1)
        throw std::runtime_error("HuffWT: file read error (Rs).");

    if (!leaf)
    {
        bitrank = new BitRank(file);
        left = new HuffWT(file, ct);
        right = new HuffWT(file, ct);
    }
}

void HuffWT::save(std::FILE *file)
{
    if (std::fwrite(&leaf, sizeof(bool), 1, file) != 1)
        throw std::runtime_error("HuffWT: file write error (Rs).");
    if (std::fwrite(&ch, sizeof(uchar), 1, file) != 1)
        throw std::runtime_error("HuffWT: file write error (Rs).");

    if (!leaf)
    {
        bitrank->save(file);
        left->save(file);
        right->save(file);
    }
}


HuffWT::~HuffWT() {
    if (left) delete left;
    if (right) delete right;
    if (bitrank)   
        delete bitrank;
}

class node 
{
private:
    unsigned weight;
    uchar value;
    node *child0;
    node *child1;
    
    void maketable( unsigned code, unsigned bits, HuffWT::TCodeEntry *codetable ) const;
    static void count_chars(uchar *, ulong , HuffWT::TCodeEntry *);
    static unsigned SetBit(unsigned , unsigned , unsigned );
public:
    node( unsigned char c = 0, unsigned i = 0 ) {
        value = c;
        weight = i;
        child0 = 0;
        child1 = 0;
    }
        
    node( node* c0, node *c1 ) {
        value = 0;
        weight = c0->weight + c1->weight;
        child0 = c0;
        child1 = c1;
    }

      
    bool operator>( const node &a ) const {
        return weight > a.weight;
    }

    static HuffWT::TCodeEntry *makecodetable(uchar *, ulong);
};




HuffWT::TCodeEntry * node::makecodetable(uchar *text, ulong n)
{
    HuffWT::TCodeEntry *result = new HuffWT::TCodeEntry[ 256 ];
    
    count_chars( text, n, result );
    std::priority_queue< node, std::vector< node >, std::greater<node> > q;
    for ( unsigned int i = 0 ; i < 256 ; i++ )
        if ( result[ i ].count )
            q.push(node( i, result[ i ].count ) );

    while ( q.size() > 1 ) {
        node *child0 = new node( q.top() );
        q.pop();
        node *child1 = new node( q.top() );
        q.pop();
        q.push( node( child0, child1 ) );
    }
    q.top().maketable(0u,0u, result);
    q.pop();
    return result;
}



void node::maketable(unsigned code, unsigned bits, HuffWT::TCodeEntry *codetable) const
{
    if ( child0 ) 
    {
        child0->maketable( SetBit(code,bits,0), bits+1, codetable );
        child1->maketable( SetBit(code,bits,1), bits+1, codetable );
        delete child0;
        delete child1;
    } 
    else 
    {
        codetable[value].code = code;    
        codetable[value].bits = bits;
    }
}

void node::count_chars(uchar *text, ulong n, HuffWT::TCodeEntry *counts )
{
    ulong i;
    for (i = 0 ; i < 256 ; i++ )
        counts[ i ].count = 0;
    for (i=0; i<n; i++)
        counts[(int)text[i]].count++; 
}

unsigned node::SetBit(unsigned x, unsigned pos, unsigned bit) {
      return x | (bit << pos);
}

HuffWT * HuffWT::makeHuffWT(uchar *bwt, ulong n)
{
    HuffWT::TCodeEntry * codetable = node::makecodetable(bwt,n);
    return new HuffWT(bwt,n, codetable, 0);
}

void HuffWT::save(HuffWT *wt,std::FILE *file)
{
    for (unsigned i = 0; i < 256; ++i)
        wt->codetable[i].save(file);
    wt->save(file);
}

HuffWT * HuffWT::load(std::FILE *file)
{
    TCodeEntry *ct = new HuffWT::TCodeEntry[ 256 ];
    for (unsigned i = 0; i < 256; ++i)
        ct[i].load(file);
    return new HuffWT(file, ct);
}

void HuffWT::deleteHuffWT(HuffWT *wt)
{
    delete [] wt->codetable;
    delete wt;
}
