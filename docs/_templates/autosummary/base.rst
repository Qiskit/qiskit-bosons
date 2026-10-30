{#
   We overwrite the standard template to enforce only showing `objname` as the
   top rather than `fullname` to avoid long page titles.
-#}

{{ objname | escape | underline}}

.. currentmodule:: {{ module }}

.. auto{{ objtype }}:: {{ objname }}
