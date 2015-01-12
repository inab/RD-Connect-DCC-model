package Heap::VCFPriority;


# This file is based on
######################################################################
#
#   Copyright (c) 2002 Frank J. Wojcik. All rights reserved.
# This program is free software; you can redistribute it and/or
#        modify it under the same terms as Perl itself.
#
######################################################################

use strict;
use Carp;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);

require Exporter;

@ISA = qw (AutoLoader);
@EXPORT = qw ();
$VERSION = "0.11";

sub new {
    
    use Heap::Fibonacci;
    use Heap::Elem::VCFElem;
    use Heap::Elem::VCFElemRev;

    my $self = shift;
    my $class = ref($self) || $self;
    my $defaults = {
	'.count'   => 0,
	'.hifirst' => 1,
	'.fifo'    => 1,
	'.queues'  => undef,
	'.heap'    => undef,     # instantiated below
	'.error'   => '',
	'.raise_error' => 0 };
    $defaults->{'.heap'} = Heap::Fibonacci->new;
    return bless $defaults, $class;
}

#
# Error handling
#

sub raise_error {
    my ($self, $errlvl) = @_;
    $errlvl ||= 0;
    unless ($errlvl =~ /^\d+$/) {
	$self->_error("Priority must be a whole number!\n");
	return;
    }
    if ( ($errlvl != 0) && ($errlvl != 1) && ($errlvl != 2) ) {
	$self->_error("Invalid error level $errlvl\n");
	return;
    }
    $self->{'.raise_error'} = $errlvl;
}

sub err_str {
    return $_[0]->{'.error'};
}

sub reset_err_str {
    return $_[0]->{'.error'} = '';
}

sub _error {
    my ($self, $error) = @_;
    $self->{'.error'} .= $error;
    croak $self->{'.error'} if $self->{'.raise_error'} == 2;
    carp  $self->{'.error'} if $self->{'.raise_error'} == 1;
}

#
# Setting of item ordering
#

sub fifo {
    $_[0]->{'.fifo'} = 1;
}

sub lifo {
    $_[0]->{'.fifo'} = 0;
}

sub lilo {
    $_[0]->{'.fifo'} = 1;
}

sub filo {
    $_[0]->{'.fifo'} = 0;
}

sub highest_first {
    my $self = $_[0];
    return if $self->{'.hifirst'} == 1;
    unless ($self->{'.count'} == 0 ) {
	my $n = Heap::Fibonacci->new;
	my $e_old;
	while(defined($e_old = $self->{'.heap'}->extract_top)) {
	    my $e_new = Heap::Elem::VCFElemRev->new;
	    $e_new->val($e_old->val);
	    $n->add($e_new);
	}
	$self->{'.heap'} = $n;	    
    }
    $self->{'.hifirst'} = 1;
}

sub lowest_first {
    my $self = $_[0];
    return if $self->{'.hifirst'} == 0;
    unless ($self->{'.count'} == 0) {
	my $n = Heap::Fibonacci->new;
	my $e_old;
	while(defined($e_old = $self->{'.heap'}->extract_top)) {
	    my $e_new = Heap::Elem::VCFElem->new;
	    $e_new->val($e_old->val);
	    $n->add($e_new);
	}
	$self->{'.heap'} = $n;	    
    }
    $self->{'.hifirst'} = 0;
}

#
# Heap information functions
#

sub count {
    return $_[0]->{'.count'};
}

sub get_priority_levels {
    my $self = shift;
    my @levels = $self->{'.hifirst'} ? sort {$b cmp $a}  keys %{$self->{'.queues'}} :
	                               sort {$a cmp $b}  keys %{$self->{'.queues'}} ;
    return wantarray ? @levels : scalar @levels;
}

sub get_level {
    my ($self, $priority) = @_;
    unless (defined($priority)) {
	$self->_error("Need to supply a priority to get_level!\n");
	return;
    }
    unless (defined($self->{'.queues'}) && exists($self->{'.queues'}->{$priority})) {
	$self->_error("Priority level $priority does not exist in heap!\n");
	return;
    }
    my @items = @{$self->{'.queues'}->{$priority}};
    @items = reverse @items unless $self->{'.fifo'};
    return wantarray ? @items : scalar @items;
}

sub get_list {
    my $self = shift;
    my @heap;
    my @levels = $self->get_priority_levels();
    push @heap, $self->get_level($_) for @levels;
    return wantarray ? @heap : scalar @heap;
}

sub get_heap {
    $_[0]->get_list();
}

sub next_level {
    return unless defined($_[0]->{'.heap'});
    my $elem = $_[0]->{'.heap'}->top;
    return $elem;
}

sub next_item {
    my $self = shift;
    my $priority = $self->next_level()->priority;
    return unless defined($priority);
    return unless defined($self->{'.queues'});
    return unless exists ($self->{'.queues'}->{$priority});
    my $idx = ($self->{'.fifo'} == 1) ? 0 : -1;
    return $self->{'.queues'}->{$priority}[$idx];
}

#
# Manipulate the heap
#

sub add {
    my ($self, $item, $priority) = @_;
    unless (defined($item)) {
	$self->_error("Need to supply an item to add to heap!\n");
	return;
    }
    $priority ||= Heap::Elem::VCFElem::genPriorityString($item);
    
    #
    # First, add the item to its priority queue
    #
    push @{$self->{'.queues'}->{$priority}}, $item;
    
    #
    # Then, add the priority to the heap if it's not already there
    #
    if (@{$self->{'.queues'}->{$priority}} == 1) {
	# For some reason, using NumRev->new($priority) DOES NOT work
	
	my $elem = ($self->{'.hifirst'} == 1) ? Heap::Elem::VCFElemRev->new :
	                                        Heap::Elem::VCFElem->new;
	$elem->val($item);
	$self->{'.heap'}->add($elem);
    }
    
    #
    # Finally, update the statistics
    #
    $self->{'.count'}++;    
}

sub pop {
    my $self = shift;
    my $priority = $self->next_level()->priority;
    return unless defined($priority);
    return unless defined($self->{'.queues'});
    my $retval = $self->{'.fifo'} ? shift @{$self->{'.queues'}->{$priority}} :
               	                    pop   @{$self->{'.queues'}->{$priority}} ;
    if (@{$self->{'.queues'}->{$priority}} == 0) {
	delete $self->{'.queues'}->{$priority};
	my $garbage = $self->{'.heap'}->extract_top;
    }
    $self->{'.count'}--;
    return $retval;
}

sub pop_queue {
	my $self = shift;
	my $priority = $self->next_level()->priority;
	return unless defined($priority);
	return unless defined($self->{'.queues'});
	
	# As we want the whole priority queue, update it
	my $retval = delete $self->{'.queues'}->{$priority};
	my $garbage = $self->{'.heap'}->extract_top;
	
	# Update the statistics
	$self->{'.count'} -= scalar(@{$retval});
	
	return $retval;
}

sub delete_item {
    my ($self, $item, $priority) = @_;
    my @levels;

    unless (defined($item)) {
	$self->_error("delete_item() must be called with an item!\n");
	return;
    }

    if (defined($priority)) {
	unless (defined($self->{'.queues'}) && exists($self->{'.queues'}->{$priority})) {
	    $self->_error("Priority $priority does not exist in heap!\n");
	    return;
	}
	@levels = ($priority);
    } 
    else {
	@levels = map { $_->priority } $self->get_priority_levels;
    }

    my $found = 0;
    for my $level (@levels) {
	my $before = @{$self->{'.queues'}->{$level}};
	@{$self->{'.queues'}->{$level}} = grep($_ ne $item, @{$self->{'.queues'}->{$level}});
	my $after =  @{$self->{'.queues'}->{$level}};

	$found += $before - $after;

	if ($after == 0) {
	    $self->delete_priority_level($level);
	}

    }

    if ($found == 0) {
	if (defined($priority)) {
	    $self->_error("Item $item does not exist at priority level $priority in heap!\n");
	} 
	else {
	    $self->_error("Item $item does not exist in heap!\n");
	}
	return;
    }
    
    $self->{'.count'} -= $found;

}

sub delete_priority_level {
    my ($self, $priority) = @_;
    unless (defined($priority)) {
	$self->_error("delete_priority_level() must be called with a priority!\n");
	return;
    }
    unless (defined($self->{'.queues'}) && exists($self->{'.queues'}->{$priority})) {
	$self->_error("Priority level $priority does not exist in heap!\n");
	return;
    }
    
    #
    # Update statistics
    #
    $self->{'.count'} -= @{$self->{'.queues'}->{$priority}};
    
    #
    # Remove that level's queue
    #
    delete $self->{'.queues'}->{$priority};
    
    #
    # Remove that level from the heap
    #
    
    # ARGH: THIS SHOULD WORK! IT DOESN'T!
    # All Heap::Elem references would need
    # to be stored in a hash or something
    # since Heap::* checks the elem references, not
    # their values! 
    #
    # Mayhaps a patch to Heap::* is in order?
    #
    #my $elem = ($self->{'.hifirst'} == 1) ? Heap::Elem::NumRev->new :
    #	                                    Heap::Elem::Num->new;
    #$elem->val($priority);
    #$self->{'.heap'}->delete($elem);

    # This is a much slower workaround. Create a new heap, and start
    # moving through priority levels (saving them on the new heap)
    # until we find the one we're deleting. Then, move all the saved
    # levels back onto the original.

    my $tmp = Heap::Fibonacci->new;
    my $i;
    while(defined($i = $self->{'.heap'}->extract_top)) {
	last if ($i->val->priority eq $priority);
	$tmp->add($i);
    }
    $self->{'.heap'}->absorb($tmp);	  
    
}

1;

__END__

=head1 NAME

Heap::Priority - Implements a priority queue or stack

=head1 SYNOPSIS

    use Heap::Priority;
    my $h = new Heap::Priority;

    $h->add($item,[$priority]); # add an item to the list
    $item = $h->pop;            # get an item back from the list

    $h->fifo;                   # set first in first out (ie. a queue) (default)
    $h->lilo;                   # synonym for fifo()
    $h->lifo;                   # set last in first out (ie. a stack)
    $h->filo;                   # synonym for lifo()
    $h->highest_first;          # set pop() in high to low priority order (default)
    $h->lowest_first;           # set pop() in low to high priority order

    $numitems  = $h->count;     # get the number of items in the list
    $priority  = $h->next_level;# get the next item's priority that would be popped
    $item      = $h->next_item; # get the next item that would be popped
    @levels    = $h->get_priority_levels;  # get the existing priority levels
    @items     = $h->get_level($priority); # get the items at a specific level
    @all_items = $h->get_list;  # get all the items on the list
    @all_items = $h->get_heap;  # synonym for get_list()

    $h->delete_item($item,[$priority]);    # delete an item from the list
    $h->delete_priority_level($priority);  # delete all items at a priority level

    $h->raise_error($err_level);# change the severity of errors
    $error_string = $h->err_str;# get the error strings
    $h->reset_err_str;          # clear the error strings

=head1 DESCRIPTION

This module implements a priority queue or stack. It will be referred
to here as a 'priority list' or 'priority heap'. The main functions
are add() and pop() which add and remove from the priority list
according to the rules you choose. When you add() an item to the
priority list you can assign a priority level to the item or let the
priority level default to 0. Multiple instances of the same item are
permitted in the priority list, either on the same or different
priority levels. A priority level must be a number, positive or
negative. Basically, any number you can fit into a scalar. If you use
floats, though, watch out for stringification issues.

What happens when you call pop() depends on the configuration you
choose. By default the highest priority items will be popped off in
First-In-First-Out order. fifo() and lifo() set First-in-First-Out and
Last-In-First-Out respectively. lilo() and filo() are synonyms for
fifo() and lifo() respectively. Calling fifo(), filo(), lifo(), or
lilo() at any time is legal, and will not change the contents of the
list.

highest_first() and lowest_first() allow you to choose to pop() the
highest priority items first or the lowest priority items first. In
the current implementation, this operation can be relatively slow,
depending on the number of different priority levels in the priority
list.

There are several functions to inspect the priority list. None of
these will change the priority list in any way. count() will return
the number of items remaining on the priority list. Information about
the next item that pop() will pull off of the priority list can be
gotten via next_level() and next_item(). get_priority_levels() will
return the list of priority levels which exist (have items assigned to
them) in the priority list. get_level() will return the list of items
at a given priority level. get_list() will return the complete set of
items in the priority list. get_heap() is a synonym for get_list();
see L<"INTERFACE COMPATIBILITY">. Items gotten via get_level(),
get_list(), or get_heap() will be returned in the order they
will/would be returned via pop().

There are two ways to manipulate the priority list, outside of the
usual pop(). delete_item() will delete all instances of an item at a
given priority level. If a priority level is not specified, then all
instances of the item at every priority level will be deleted.  This
can be very slow. See L<"EFFICIENCY">. delete_priority_level() will
delete all items at a given priority level. These methods are not
deprecated, but work "outside" of the usual purpose of a priority
list. If you find your program using them, perhaps you should consider
using a different data structure.

Error handling is modified via raise_error(). When an error is
generated, the default is to record it, and continue on silently. This
can be changed by calling raise_error() with an error level. If the
error level is set to 1, errors will be warned about (via carp()). If
the error level is set to 2, errors will be fatal (via croak()).
Default behavior can be reset by setting the error level to 0. Error
messages produced can be accessed through err_str(). Error messages
can be cleared by calling reset_err_str().

=head1 INTERFACE COMPATIBILITY

This interface is designed to be almost 100% drop-in compatible with
Heap::Priority 0.01 by Dr James Freeman E<lt>jfreeman@tassie.net.auE<gt>. 
About 20-30% of the code and documentation is extremely similar to
the original package.

The major exception is that the modify_priority() function is not
implemented. It is not currently possible to change the priority of an
item. To workaround this, you can explicitly do what the previous
implementation implicitly did, which is to delete all instances of the
item and add in a new one (via delete_item() and add()).

Further, the terminology of that package was somewhat looser than I
would have liked, hence the synonyms of 'list' and 'heap'. The error
messages returned by err_str() are almost identical; one "on heap!"
was changed to an "in heap!". Exact error texts should not be relied
upon to stay the same in future revisions.

=head1 OBJECT INTERFACE

This is an OO module. You begin by creating a new heap object

    use Heap::Priority;
    my $h = new Heap::Priority;

You then simply call methods on your heap object:

    $h->lowest_first();         # set pop() in low to high priority order
    $h->add($item,$priority);   # add an item to the heap
    $h->lifo;                   # set last in first out (ie. a stack)
    my $item = $h->pop;         # get an item back from heap

=head1 METHODS

=head2 new()

    my $h = new Heap::Priority;

The constructor takes no arguments and simply returns an empty default
list. The default configuration is FIFO (ie. a queue) with highest
priority items popped first

=head2 add($item,[$priority])

    $h->add($item, [$priority]);

add() will add $item to the heap. Optionally, it may be assigned a
$priority level (default priority level is 0).

=head2 pop()

    my $item = $h->pop;

pop() takes no arguments. In default configuration pop() will return
those items having the highest priority level first in FIFO
order. This behavior can be modified using the methods outlined below.

=head2 fifo()

    $h->fifo;

Set pop() to work on a First-In-First-Out basis, otherwise known as a
queue.  This is the default configuration.

=head2 lilo()

    $h->lilo;

This is a synonym for fifo().

=head2 lifo()

    $h->lifo;

Set pop() to work on a Last-In-First-Out basis, otherwise known as a
stack.

=head2 filo()

    $h->filo;

This is a synonym for lifo().

=head2 highest_first()

    $h->highest_first;

Set pop() to retrieve items in highest to lowest priority order. This
is the default configuration. 

=head2 lowest_first()

    $h->lowest_first;

Set pop() to retrieve items in lowest to highest priority order. 

=head2 count()

    my $numitems = $h->count;

Returns the number of items in the priority list. Much lighter weight
than counting the items returned by get_list().

=head2 next_level()

    my $priority = $h->next_level;

Returns the priority level assigned to the next item that pop() would
return. Does not modify the priority list.

=head2 next_item()

    my $item = $h->next_item;

Returns the next item that pop() would return. Does not modify the
priority list. (Otherwise, it would be identical to pop())

=head2 get_priority_levels()

    my @levels = $h->get_priority_levels;

Returns the list of priority levels in current pop() order.

=head2 get_level($priority)

    my @items = $h->get_level($priority);

Returns the entire list of items assigned to the specified priority
level in current pop() order.

=head2 get_list()

    my @all_items = $h->get_list;

Returns entire list of items in current pop() order. If you only want
the number of items on the list, count() is probably faster.

=head2 get_heap()

    my @all_items = $h->get_heap;

Synonym for get_list(). See L<"INTERFACE COMPATIBILITY">. 

=head2 delete_item($item,[$priority])

    $h->delete_item($item,[$priority]);

This method will delete $item from the heap. If the optional $priority
is not supplied all instances of $item will be removed from the
heap. If $priority is supplied then only instances of $item at that
priority level will be removed.

=head2 delete_priority_level($priority)

    $h->delete_priority_level($priority);

Delete all items assigned to $priority level.

=head2 raise_error($n)

    $h->raise_error(1);

Set error level $n => 2 = croak, 1 = carp, 0 = silent (default)

=head2 err_str()

    my $error_string = $h->err_str;

Return error string, if any.

=head2 reset_err_str()

    $h->reset_err_str;

Reset any errors that err_str() would have reported.

=head1 EXPORT

Nothing: it's an OO module.

=head1 IMPLEMENTATION DETAILS

Items are stored in lists, one list per priority level. These lists
are stored in a hash, with the priority level as the key. The priority
levels are stored in a true heap (via L<Heap::Fibonacci>,
L<Heap::Num>, and L<Heap::NumRev>). Therefore, theoretically,
arbitrary scalars could be used as priority levels with a minimal
amount of code change.

Currently, only one B<Heap::Fibonacci> instance is stored
persistently. The behavior when highest_first and lowest_first are
called on non-empty priority lists is to create a new instance, and
move the elements to the new heap, converting their type as
appropriate (B<Heap::Num> to/from B<Heap::NumRev>). It would also be
possible to concurrently store two heaps, one for highest_first and
one for lowest_first. I think I've decided that the first option is
better, because changing that particular horse in mid-stream is
uncommon enough to justify the memory savings in the more usual
cases. I could be convinced this is not true. As it stands, you should
probably decide on priority list behavior before adding items.

The reason that modify_priority() is not implemented is that an item's
position relative to the others' is not recorded, so strict fifo/lifo
behavior could not be guaranteed. This is completely aside from the
problem of selecting and specifying which instance of the item to
modify, though "all of them" seems a not unreasonable behavior.

There is no test for item existence, no way to count the instances of
a particular item on the list, no way to list the priorities an item
is assigned to, no way delete specific instances of an item that is on
the list multiple times, and no way to query the priority list about
its current pop() ordering. Implementing these would require some
combination of changing the interface, using a much more memory-hungry
implementation, or making those calls very, very slow. I am
considering writing a heavier-weight (more memory intensive) version
of this called B<Heap::Priority::Super> (or some-such) which does
support those features. Another way to solve this problem would be to
pass an optimization parameter on object instantiation. This may be
simpler from a user point of view, but I think the code would become a
mess of C<if> statements. Comments are very welcome.

=head1 EFFICIENCY

The internal data storage representation was chosen to maximize speed
for "usual" operations, while minimizing memory usage. List
modification outside of add() and pop() might be slow; in particular,
deleting an item completely from the list will require visiting every
item on the list.

The representation is efficient for very large numbers of items,
regardless of the number of distinct priority levels.

Benchmarks (on a 750 MHz Pentium III, running FreeBSD 4.5 and perl
5.00503) indicate that add() and pop() take about 25 usec,
best-case. Real-world testing on the same machine, with approximately
1,000,000 items 10,000 priority levels shows add() taking about 55
usec, and pop() taking about 50 usec, on average.  It turns out that
perl subroutine/method overhead is substantial compared to the rest of
the program, and I was able to obtain significant speedups by
'reaching in' to the objects, bypassing accessor methods. This is
B<not> done in this release of the code. A future implementation may
not use B<Heap::Fibonacci>, but may hard-code a number-based heap
implementation internally.

A complexity analysis of this code, assuming n items evenly distributed
over m levels:

O(1)

=over 4

=item new()       

=item fifo()/lilo()

=item filo()/lifo()

=item count()

=item next_level() 

=item next_item()  

=back

O(log m)

=over 4

=item add()

=item pop()

=back

O(m/n)

=over 4

=item get_level()

=back

O(m log m)

=over 4

=item highest_first()

=item lowest_first()

=item get_priority_levels()

=item delete_priority_level()

=back

O(m log m + n/m)

=over 4

=item delete_item(), specific level

=back

O(m log m + n)

=over 4

=item get_list()/get_heap()

=back

O(m^2 log m + n)

=over 4

=item delete_item(), all levels         

=back

=head1 TODO

Add better tests for get_priority_levels() and get_level(). Add more
tests with multiple items on the same level to test fifo/lifo
correctness.

=head1 BUGS

There's always one bug. If you find it, please let me know.

=head1 SEE ALSO

Heap(3), Heap::Fibonacci(3), Heap::Elem::Num(3), Heap::Elem::NumRev(3)

=head1 AUTHOR

Frank J Wojcik E<lt>fwojcik+pri@besh.comE<gt>

=cut
